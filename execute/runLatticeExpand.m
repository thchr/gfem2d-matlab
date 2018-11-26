function savepath = runLatticeExpand(latticetype,kspec,Ns,Ncirc,adivd,expandshrink,Nk,savephi,models,addstr)

fprintf('\n2D lattice calculation commenced on | %s\n',datestr(now));
addpath(genpath('..'))

%%  ALLOW MULTITHREADING ON CLUSTER
fprintf('\n|----- MULTITHREADING SETUP -----|\n'); nopool = 1; %Don't open a multiworker pool
SetMatlabMultithreading;

%% SETUP DEFAULTS
if nargin < 1 || isempty(latticetype)
    latticetype = 'triangular'; %Lattice type (square is the alternative)
end
if nargin < 2 || ~exist('kspec','var') || isempty(kspec)
    kspec = 'irrfbz'; %Number of points along each unit cell side
end
if nargin < 3 || ~exist('Ns','var') || isempty(Ns)
    Ns = 118; %Number of points along each unit cell side
end
if nargin < 4 || ~exist('Ncirc','var') || isempty(Ncirc)
    Ncirc = 70; %Number of points used in circular inclusion
end
if nargin < 5 || ~exist('adivd','var') || isempty(adivd)
    adivd = 2; %Period divided by inclusion diameter
end
if nargin < 6 || ~exist('expandshrink','var') || isempty(expandshrink)
    expandshrink = 0; %factor to expand (+) or shrink (-) the hexagonal "atom" with
end
if nargin < 7 || ~exist('Nk','var') || isempty(Nk)
    if ischar(kspec);  switch kspec                %Number of k-points
            case 'irrfbz'; Nk = 101;  case 'projx';  Nk = [49,11];
        end
    else  Nk = []; end
end
if nargin < 8 || ~exist('savephi','var') || isempty(savephi)
    savephi = 0; %Whether to save the potential phi (1) or not (0)
end
if nargin < 10 || ~exist('addstr','var') || isempty(addstr)
    addstr = []; %Whether to add a modifier string to the savename ([] = nothing added)
elseif ~isempty(addstr)
    addstr = ['_' addstr];
end 

%% PRINT SETUP
if iscellstr(kspec) %Expect { kx, ky } with kx,ky entries as _strings_
    if numel(kspec) == 2;  kspec{3} = 'gridded';  end %Default; gridded (otherwise; list)
    %If kspec is supplied as a cell, we need to take this
    %into account for the following fprint command
    kspecstr = cellstr(kspec);
    kspecstr = cellfun(@(x) cat(2,x,' & '),kspecstr,'UniformOutput',false);
    kspecstr = [kspecstr{:}]; kspecstr = kspecstr(1:end-3);
else
    kspecstr = kspec; %If kspec is not a cell array
end

strf = @(instr,strtype) [repmat(' ',1,14-numel(instr)) instr ' | %' strtype '\n']; %Convenient print string as function

fprintf('\n|----- SETUP -----|\n')
fprintf([strf('latticetype','s'), strf('kspec','s'), strf('Ns','g'), strf('Ncirc','g'), strf('adivd','g'), strf('expandshrink','g'), strf('Nk','s')              , strf('savephi','g'), strf('addstr','s')],...
         latticetype            , kspecstr         , Ns            , Ncirc            , adivd            , expandshrink,             strrep(num2str(Nk),'  ','x'), savephi            , addstr)

%% MESHING & MOMENTA
fprintf('\n|----- MESHING UNIT CELL -----|\n')
theta = linspace(0,2*pi,Ncirc+1); theta = theta(1:end-1); %Inclusion boundary
hole = 1/(2*adivd)*[cos(theta);sin(theta)].';        %----||----

R = calcDirect(latticetype);

inclusions=cell(1,4); 
cntrs = [-1,-1; 1,-1; 1,1; -1,1]/4; cntrs = cntrs(:,1)*R{1} + cntrs(:,2)*R{2};
[cntrstheta,cntrsrho] = cart2pol(cntrs(:,1),cntrs(:,2));
newcntrs = bsxfun(@times, [cos(cntrstheta),sin(cntrstheta)], cntrsrho*(1+expandshrink));
for ii = 1:4
    inclusions{ii} = bsxfun(@plus, hole/2, newcntrs(ii,:));
end

mesh = geomPeriodicInclusion(latticetype,inclusions,Ns,@(x,y) .015*250/Ns*ones(size(x)),1);

%Momentum list
if ischar(kspec)
    switch kspec
        case 'projx'
            [kx,ky] = meshgrid(linspace(0,pi,Nk(1)),linspace(0,1,Nk(2))*2*pi/sqrt(3));
            kvec.x = kx.'; kvec.y = ky.'; Nk = Nk(1)*Nk(2);        clear kx ky
        case 'irrfbz'
            [kxy,kplot,kmark] = irreducibleFBZ('triangular',Nk(1),0);
            kvec.x = kxy(:,1); kvec.y = kxy(:,2);                  clear kxy
    end
elseif iscellstr(kspec)
    
    kx = eval(kspec{1}); ky = eval(kspec{2});
    if strcmpi(kspec{3},'gridded')  %In this case, we interpret kspec to be intended as a gridded list
        [kx,ky] = meshgrid(kx,ky);
    end
    kvec.x = kx.'; kvec.y = ky.'; Nk = numel(kx);          clear kx ky
    kspecid = kspecstr;
    kspecstr = 'manualkspec';
end

fprintf('|----- REQUESTED MOMENTA (%g k-points) -----|\n',Nk)
switch kspecstr
    case 'manualkspec'
        fprintf('Manual momentum-grid requested, according to command:\n')
        fprintf('   %s\n   [kx] = [',kspecid)
        fprintf('%-.3f,',unique(kvec.x)); fprintf('\b]\n   [ky] = ['); fprintf('%-.3f,',unique(kvec.y)); fprintf('\b]');
    case 'projx'
        fprintf('Square momentum-grid spanned by k-point [kx] x [ky] grid:\n   [kx] = [')
        fprintf('%-.3f,',unique(kvec.x)); fprintf('\b]\n   [ky] = ['); fprintf('%-.3f,',unique(kvec.y)); fprintf('\b]');
    case 'irrfbz'
        fprintf('First BZ spanned by k-point pairs:\n   [kx,ky] = ')
        fprintf('[%-.3f,%-.3f] ',[kvec.x(:),kvec.y(:)].')
end


%% MODELS
if nargin < 9 || ~exist('models','var') || isempty(models) %Standard, default model
%% SET UP A NUMBER OF MODELS
    ef_eV = .4; L = 400e-9;
    clear models %Apparently, it is necessary (???) to remove variable 'models' from namespace, even if it is empty, in order to assign it as a struct with entries
    models(1) =     struct('type','isotropic','ef_eV',ef_eV,'L',L,'approx','intra','keep_eV',[0,ef_eV]);
    %models(end+1) = struct('type','isotropic','ef_eV',ef_eV,'L',L,'approx','full' ,'keep_eV',[0,ef_eV]);
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
fprintf('\n\n%g ''setup-models'' are considered:\n',numel(models))
modelsprint = struct2table(models); disp(modelsprint); clear modelsprint

%% ACTUAL CALCULATION AND PLOTTING
bloch.type = '2D';

fprintf('\n\n|----- LATTICE DISPERSION CALCULATION -----|\n')
tick = tic;
for kk2 = size(kvec.x,2):-1:1 %This loop is only relevant in case kspec = 'irrfbz'
    for kk1 = size(kvec.x,1):-1:1
        %--- ACTUAL CALCULATION ---
        bloch.k = [kvec.x(kk1,kk2),kvec.y(kk1,kk2)];
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
            eneeig_eV(1:numel(eneeig_eV_cell),kk1,kk2) = eneeig_eV_cell;
            if savephi >= 1; phieig{kk1,kk2} = phi_cell; end
            if savephi == 2; rhoeig{kk1,kk2} = rho_cell; end
        else %Otherwise, we group in cell so that each model is in its own cell
            for mm = 1:numel(models)
                eneeig_eV{mm}(1:numel(eneeig_eV_cell{mm}),kk1,kk2) = eneeig_eV_cell{mm};
                if savephi >= 1; phieig{mm}{kk1,kk2} = phi_cell{mm}; end
                if savephi == 2; rhoeig{mm}{kk1,kk2} = rho_cell{mm}; end
            end
        end
        
        %--- PRINT PROGRESS ---
        fprintf('--- DISPERSION LOOP %g/%g (%.2f/%.2f min) ---\n\n',Nk+1-kk1-(kk2-1)*size(kvec.x,1),Nk,toc(tick)/60,toc(tick)/60*Nk/(Nk+1-kk1-(kk2-1)*size(kvec.x,1)));
    end
end

%Now we want to assign all "high" (row-sense) elements in eneeig_eV{mm}
%to NaN if they are currently zero (i.e. they are unassigned values and not
%a valid zero mode; i.e. they do not reside at kvec = [0,0])
for mm = 1:numel(models)
    for kk2 = size(kvec.x,2):-1:1
        %Find the elements of eneeig_eV which equal zero
        if numel(models) == 1
            [rowz,colz] = find(eneeig_eV(:,:,kk2) == 0);
        else
            [rowz,colz] = find(eneeig_eV{mm}(:,:,kk2) == 0);
        end
        %Set any =0-element to NaN, unless it is the first band at k=[0,0]
        for bb = rowz(:).'          %Matlab iterates over columns, so make sure this is a row vector
            if bb ~= 1
                for kk1 = colz(:).' %Matlab iterates over columns, so make sure this is a row vector
                    if kvec.x(kk1,kk2) ~= 0 || kvec.y(kk1,kk2) ~= 0
                        if numel(models) == 1
                            eneeig_eV(bb,kk1,kk2) = NaN;
                        else
                            eneeig_eV{mm}(bb,kk1,kk2) = NaN;
                        end
                    end
                end
            end
        end
    end
end


%% SAVE DATA
savedir = '../output/lattice/';
savename = ['lattice_expandedtriangular_' strrep(strrep(num2str(expandshrink*100,'%+g'),'-','m'),'+','p') '_adivd' num2str(adivd) '_' kspecstr '_Ns' num2str(Ns) '_Ncirc' num2str(Ncirc) '_savephi' num2str(savephi) addstr];
savepath = [savedir savename '.mat'];
clear eneeig_eV_cell phi_cell rho_cell tick mm kk strf poolobj numProc_cluster numProc_private LASTppn CurMatPool savename savedir mm kk1 kk2 bb

save(savepath)

%Print final messages
fprintf('\n\n All data saved to   | %s\n\n', savepath)
fprintf('Calculated finished on   | %s\n  ', datestr(now))