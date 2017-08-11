function savepath = runRibbonLatticeMagneto(latticetype,numcell,Ns,Ncirc,adivd,Nk,savephi,models)

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
if nargin < 3 || ~exist('Ns','var') || isempty(Ns)
    Ns = 70; %Number of points along each unit cell side
end
if nargin < 4 || ~exist('Ncirc','var') || isempty(Ncirc) 
    Ncirc = 30; %Number of points used in circular inclusion
end 
if nargin < 5 || ~exist('adivd','var') || isempty(adivd) 
    adivd = 2; %Period divided by inclusion diameter
end  
if nargin < 6 || ~exist('Nk','var') || isempty(Nk)
    Nk = 49; %Number of k-points
end 
if nargin < 7 || ~exist('savephi','var') || isempty(savephi)
    savephi = 0; %Whether to save the potential phi (1) or not (0)
end 

%Print setup
fprintf('\n|----- SETUP -----|\n')
strf = @(instr,strtype) [repmat(' ',1,14-numel(instr)) instr ' | %' strtype '\n']; %Convenient print string as function
fprintf([strf('lattice type','s'), strf('numcell','g'), strf('Ns','g'), strf('Ncirc','g'), strf('adivd','g'), strf('Nk','g'), strf('savephi','g')],...
         latticetype             ,numcell             , Ns            , Ncirc           , adivd             , Nk            , savephi)

%% MESHING & MOMENTA
fprintf('\n|----- MESHING -----|\n')
theta = linspace(0,2*pi,Ncirc+1); theta = theta(1:end-1); %Inclusion boundary
inclusion = 1/(2*adivd)*[cos(theta);sin(theta)].';        %----||----
mesh = geomRibbonInclusion(latticetype, [1,numcell],inclusion,Ns,...
                           @(x,y) .015*250/Ns*ones(size(x)),0);

%Momentum list
kvec = bsxfun(@times,linspace(-1/2,1/2,Nk)',mesh.Grib); 
maxk = max(sqrt(kvec(:,1).^2 + kvec(:,2).^2));

%% MODELS 
ef_eV = .2; L = 100e-9; 
if nargin < 8 || ~exist('models','var') || isempty(models) %Standard, default model
    models(1) = struct('type','isotropic','ef_eV',ef_eV,'L',L,'B',[],'system',[]        ,'approx','intra','keep_eV',[0,ef_eV]);
    models(2) = struct('type','magneto'  ,'ef_eV',ef_eV,'L',L,'B',5 ,'system','graphene','approx',[]     ,'keep_eV',[.0135,ef_eV]); %The lower cutoff is the cyclotron frequency
    models(3) = struct('type','magneto'  ,'ef_eV',ef_eV,'L',L,'B',10,'system','graphene','approx',[]     ,'keep_eV',[.0270,ef_eV]); %---||---
end
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
        [eneeig_eV_cell,~,phi_cell] = calcEigenEnergiesAny(mesh,models,bloch);
    end
    for mm = 1:numel(models)
        eneeig_eV{mm}(1:numel(eneeig_eV_cell{mm}),kk) = eneeig_eV_cell{mm}; 
        if savephi == 1 
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
clear eneeig_eV_cell phi_cell tick mm kk strf poolobj numProc_cluster numProc_private LASTppn CurMatPool
savedir = '../output/ribbons/';
savename = ['magnet_' latticetype '_numcell' num2str(numcell) '_Ns' num2str(Ns) '_Ncirc' num2str(Ncirc) '_savephi' num2str(savephi)];
savepath = [savedir savename '.mat'];
save(savepath)

%Print final messages
fprintf('\n\n All data saved to   | %s\n\n', savepath)
fprintf('Calculated finished on   | %s\n  ', datestr(now))