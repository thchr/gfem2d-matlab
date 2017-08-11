function savepath = runRibbonScissorLatticeMagneto(latticetype,numcell,cut,Ns,Ncirc,adivd,Nk,savephi,models,addstr)

fprintf('\nCalculation commenced on | %s\n',datestr(now)); 
addpath(genpath('..'))

%%  ALLOW MULTITHREADING ON CLUSTER
fprintf('\n|----- MULTITHREADING SETUP -----|\n')
SetMatlabMultithreading;

%% SETUP DEFAULTS
if nargin < 1 || isempty(latticetype)
    latticetype = 'triangular'; %Lattice type (square is the alternative)
end
if nargin < 2 || ~exist('numcell','var') || isempty(numcell)
    numcell = 3; %Number of repeated unit cells in ribbon (vertical dir.)
end
if nargin < 3 || ~exist('cut','var') || isempty(cut)
	cut = 1/2; %Reduction of edge unit cells in percent along R2 
end
if nargin < 4 || ~exist('Ns','var') || isempty(Ns)
    Ns = 70; %Number of points along each unit cell side
end
if nargin < 5 || ~exist('Ncirc','var') || isempty(Ncirc) 
    Ncirc = 36; %Number of points used in circular inclusion
end 
if nargin < 6 || ~exist('adivd','var') || isempty(adivd) 
    adivd = 2; %Period divided by inclusion diameter
end  
if nargin < 7 || ~exist('Nk','var') || isempty(Nk)
    Nk = 49; %Number of k-points
end 
if nargin < 8 || ~exist('savephi','var') || isempty(savephi)
    savephi = 0; %Whether to save the potential phi (1) or not (0)
end 
if nargin < 10 || ~exist('addstr','var') 
    addstr = []; %Whether to add a modifier string to the savename ([] = nothing added)
elseif ~isempty(addstr)
    addstr = ['_' addstr];
end 

%Print setup
fprintf('\n|----- SETUP -----|\n')
strf = @(instr,strtype) [repmat(' ',1,14-numel(instr)) instr ' | %' strtype '\n']; %Convenient print string as function
fprintf([strf('lattice type','s'), strf('numcell','g'), strf('cut','g') , strf('Ns','g'), strf('Ncirc','g'), strf('adivd','g'), strf('Nk','s')          , strf('savephi','g'), strf('addstr','s')],...
         latticetype             , numcell            , cut             , Ns            , Ncirc            , adivd            , num2str(Nk(:).','%.3f,'), savephi            , addstr)

%% MESHING & MOMENTA
fprintf('\n|----- MESHING -----|\n')
theta = linspace(0,2*pi,Ncirc+1); theta = theta(1:end-1); %Inclusion boundary
inclusion = 1/(2*adivd)*[cos(theta);sin(theta)].';        %----||----
mesh = geomRibbonInclusionScissor(latticetype, numcell,inclusion,cut,Ns,...
                           @(x,y) .015*250/Ns*ones(size(x)),0);

%Momentum list
if numel(Nk) == 1
    kvec = bsxfun(@times,linspace(-1/2,1/2,Nk)',mesh.Grib);
elseif numel(Nk) > 1 %Interpret Nk as a list of points in the range [-1,1], being multiplied onto mesh.Grib
    kvec = bsxfun(@times,Nk(:),mesh.Grib); Nk = size(kvec,1);
end

%% MODELS 
if nargin < 9 || ~exist('models','var') || isempty(models) %Standard, default model
    ef_eV = .2; L = 400e-9; omegac_eV = @(B) 5.403927588936477e-04*B/ef_eV;
    clear models %Apparently, it is necessary (???) to remove variable 'models' from namespace, even if it is empty, in order to assign it as a struct with entries
    models(1) = struct('type','isotropic','ef_eV',ef_eV,'L',L,'B',0,'system',[]        ,'approx','intra','keep_eV',[0,.6*ef_eV]);
    models(2) = struct('type','magneto'  ,'ef_eV',ef_eV,'L',L,'B',4,'system','graphene','approx',[]     ,'keep_eV',[0,.6*ef_eV]); %The lower cutoff is the cyclotron frequency
    models(3) = struct('type','magneto'  ,'ef_eV',ef_eV,'L',L,'B',8,'system','graphene','approx',[]     ,'keep_eV',[0,.6*ef_eV]); %---||---
elseif ischar(models) %Evaluate a script which describes a set of models
    modelscall = models; clear models; 
    if any(modelscall == '(') %Crude way of checking if string is a function call
        models = eval(modelscall);
    else %Otherwise, we assume the script is a script, which will output models itself
        eval(modelscall); 
    end
    clear modelscall
end
%Print model type
fprintf('The following model-parameters were input:\n')
modelsprint = struct2table(models); disp(modelsprint); clear modelsprint
%% ACTUAL CALCULATION AND PLOTTING
bloch.type = '1D'; bloch.n_genexpint = 10; 

fprintf('|----- DISPERSION CALCULATION (%g k-points) -----|\n',Nk)
tick = tic;
for kk = Nk:-1:1
    %--- ACTUAL CALCULATION ---
    bloch.k = kvec(kk,:);
    if savephi == 0 %Only save the eigenenergies
        eneeig_eV_cell = calcEigenEnergiesAny(mesh,models,bloch);
    elseif savephi == 1 %Also save the eigenpotentials
        [eneeig_eV_cell,rho_cell,phi_cell] = calcEigenEnergiesAny(mesh,models,bloch);
    end
    for mm = 1:numel(models)
        eneeig_eV{mm}(1:numel(eneeig_eV_cell{mm}),kk) = eneeig_eV_cell{mm}; 
        if savephi == 1 
            rhoeig{mm}{kk} = rho_cell{mm};
            phieig{mm}{kk} = phi_cell{mm};
        end
    end
    %--- PRINT PROGRESS ---
    fprintf('--- DISPERSION LOOP %g/%g (%.2f/%.2f min) ---\n\n',Nk+1-kk,Nk,toc(tick)/60,toc(tick)/60*Nk/(Nk+1-kk));
end

%Now we want to assign all "high" (row-sense) elements in eneeig_eV{mm}
%to NaN if they are currently zero (i.e. they are unassigned values and not
%a valid zero mode)
for mm = 1:numel(models)
    [rowz,colz] = find(eneeig_eV{mm} == 0);
    eneeig_eV{mm}(sub2ind(size(eneeig_eV{mm}),rowz(rowz~=1),colz(rowz~=1))) = NaN; %This sets any element whose row is larger than 1 and whose value is equal to zero to NaN
end

%% SAVE DATA

savedir = '../output/ribbons/scissored/';
savename = ['scmag_' latticetype '_numcell' num2str(numcell) '_cut' num2str(cut*100,2) '_Ns' num2str(Ns) '_Ncirc' num2str(Ncirc) '_savephi' num2str(savephi) addstr];
savepath = [savedir savename '.mat'];
clear eneeig_eV_cell phi_cell rho_cell tick mm kk strf poolobj numProc_cluster numProc_private LASTppn CurMatPool savename savephi savedir rowz colz 
save(savepath)

%Print final messages
fprintf('\n\n All data saved to   | %s\n\n', savepath)
fprintf('Calculated finished on   | %s\n  ', datestr(now))