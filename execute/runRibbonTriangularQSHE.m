function savepath = runRibbonTriangularQSHE(adivd12,numcell,cut,Ns,Ncirc,Nk,savephi,models,addstr)
clc
fprintf('\nCalculation commenced on | %s\n',datestr(now));
addpath(genpath('..'))

%%  ALLOW MULTITHREADING ON CLUSTER
fprintf('\n|----- MULTITHREADING SETUP -----|\n')
SetMatlabMultithreading;

%% SETUP DEFAULTS
if nargin < 1 || isempty(adivd12)
    adivd12 = [1.9,2.4]; %Period divided by inclusion diameter
end
if nargin < 2 || ~exist('numcell','var') || isempty(numcell)
    numcell = 30; %Number of repeated unit cells in ribbon (vertical dir.)
end
if nargin < 3 || ~exist('cut','var') || isempty(cut)
    cut = 0; %Reduction of edge unit cells in percent along R2
end
if nargin < 4 || ~exist('Ns','var') || isempty(Ns)
    Ns = 90; %Number of points along each unit cell side
end
if nargin < 5 || ~exist('Ncirc','var') || isempty(Ncirc)
    Ncirc = [71,56]; %Number of points used in circular inclusion
end
if numel(Ncirc) == 1; Ncirc = repmat(Ncirc,size(adivd12)); end %duplicate if necessary (to match adiv12)
if nargin < 6 || ~exist('Nk','var') || isempty(Nk)
    Nk = 49; %Number of k-points
end
if nargin < 7 || ~exist('savephi','var') || isempty(savephi)
    savephi = 2; %Whether to save the potential phi (1) or not (0)
end
if nargin < 9 || ~exist('addstr','var')
    addstr = []; %Whether to add a modifier string to the savename ([] = nothing added)
elseif ~isempty(addstr)
    addstr = ['_' addstr];
end

%% LOAD MESH (?)
%Allow user to manually specify a meshfile in adiv12 (which then overrides 
%adiv12; specified relative to runRibbonTriangularQSHE). Loaded data should
%also specify adivd12, Ns, Ncirc (as array), numcell, and cut (in addition
%to a valid mesh struct (created e.g. by geomRibbonDistinctInclusions(..)]
if ischar(adivd12) 
    temploadname = adivd12; clear adivd12
    load([temploadname '.mat'])
    loadedMesh = 1;
    fprintf('\n|----- MESH LOADED FROM: %s -----|\n',temploadname); clear temploadname
else
    loadedMesh = 0;
end

%% PRINT SETUP
fprintf('\n|----- SETUP -----|\n')
strf = @(instr,strtype) [repmat(' ',1,14-numel(instr)) instr ' | %' strtype '\n']; %Convenient print string as function
fprintf([strf('lattice type','s'), strf('numcell','g'), strf('cut','g') , strf('Ns','g'), strf('Ncirc','s')        , strf('adivd','s') '\b\b\n'   , strf('Nk','s') '\b\b\n' , strf('savephi','g'), strf('addstr','s')],...
         'triangular'            , numcell            , cut             , Ns            , num2str(Ncirc(:).','%d,'), num2str(adivd12(:).','%.3f,'), num2str(Nk(:).','%.3f,'), savephi            , addstr)

%% MESHING
if ~loadedMesh
    fprintf('\n|----- MESHING -----|\n')
    
    for a=1:numel(adivd12)
        theta = linspace(0,2*pi,Ncirc(a)+1); theta = theta(1:end-1); %Inclusion boundary
        inclusion{a} = 1/(2*adivd12(a))*[cos(theta);sin(theta)].';        %----||----
    end
    mesh = geomRibbonDistinctInclusions('triangular', numcell,inclusion,cut,Ns,...
        @(x,y) .045*ones(size(x)),0);
end

%% MOMENTA
%Momentum list
if numel(Nk) == 1
    kvec = bsxfun(@times,linspace(-1/2,1/2,Nk)',mesh.Grib);
elseif numel(Nk) > 1 %Interpret Nk as a list of points in the range [-1,1], being multiplied onto mesh.Grib
    kvec = bsxfun(@times,Nk(:),mesh.Grib); Nk = size(kvec,1);
end

%% MODELS
if nargin < 8 || ~exist('models','var') || isempty(models) %Standard, default model
    ef_eV = .4; L = 400e-9;
    clear models %Apparently, it is necessary (???) to remove variable 'models' from namespace, even if it is empty, in order to assign it as a struct with entries
    models(1) = struct('type','isotropic','ef_eV',ef_eV,'L',L,'approx','intra','keep_eV',[0,.5*ef_eV]);
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
        [eneeig_eV_cell,~,phi_cell] = calcEigenEnergiesAny(mesh,models,bloch);
    elseif savephi == 2 %Save both eigenpotentials and eigendensitites
        [eneeig_eV_cell,rho_cell,phi_cell] = calcEigenEnergiesAny(mesh,models,bloch);
    end
    
    %If there's only one model, there's no need to have an extra cell-
    %index (and the output is also not a cell, just an array)
    if numel(models)==1
        eneeig_eV(1:numel(eneeig_eV_cell),kk) = eneeig_eV_cell;
        if savephi >= 1; phieig{kk} = phi_cell; end
        if savephi == 2; rhoeig{kk} = rho_cell; end
    else %Otherwise, we group in cell so that each model is in its own cell
        for mm = 1:numel(models)
            eneeig_eV{mm}(1:numel(eneeig_eV_cell{mm}),kk) = eneeig_eV_cell{mm};
            if savephi >= 1; phieig{mm}{kk} = phi_cell{mm}; end
            if savephi == 2; rhoeig{mm}{kk} = rho_cell{mm}; end
        end
    end
    
    %--- PRINT PROGRESS ---
    fprintf('--- DISPERSION LOOP %g/%g (%.2f/%.2f min) ---\n\n',Nk+1-kk,Nk,toc(tick)/60,toc(tick)/60*Nk/(Nk+1-kk));
end

%Now we want to assign all "high" (row-sense) elements in eneeig_eV{mm}
%to NaN if they are currently zero (i.e. they are unassigned values and not
%a valid zero mode; i.e. they do not reside at kvec = [0,0])
for mm = 1:numel(models)
    for kk = Nk:-1:1
        %Find the elements of eneeig_eV which equal zero
        if numel(models) == 1
            rowz = find(eneeig_eV(:,kk) == 0);
        else
            rowz = find(eneeig_eV{mm}(:,kk) == 0);
        end
        %Set any =0-element to NaN, unless it is the first band at k=[0,0]
        for bb = rowz(:).' %Matlab iterates over columns, so make sure this is a row vector
            if bb ~= 1 && (kvec(kk,1) ~= 0 || kvec(kk,2) ~= 0)
                if numel(models) == 1
                    eneeig_eV(bb,kk) = NaN;
                else
                    eneeig_eV{mm}(bb,kk) = NaN;
                end
            end
        end
    end
end

%% SAVE DATA

savedir = '../output/ribbons/qshe/';
savename = ['ribbon_qshe_triangular_adivd12_' strrep(num2str(adivd12),'         ','-') ...
            'numcell' num2str(numcell) '_Ns' num2str(Ns) '_Ncirc' strrep(num2str(Ncirc),'  ','-') ...
            '_savephi' num2str(savephi) addstr];
        
savepath = [savedir savename '.mat'];
clear eneeig_eV_cell phi_cell rho_cell tick mm kk strf poolobj numProc_cluster numProc_private LASTppn CurMatPool savename savephi savedir rowz theta inclusions
save(savepath)

%Print final messages
fprintf('\n\n All data saved to   | %s\n\n', savepath)
fprintf('Calculated finished on   | %s\n  ', datestr(now))