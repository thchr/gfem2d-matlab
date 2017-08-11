function savepath = runDiskLatticeMagneto(latticetype,kspec,Ns,Ncirc,adivd,Nk,savephi,models,savename)

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
    Ns = 70; %Number of points along each unit cell side
end 
if nargin < 4 || ~exist('Ncirc','var') || isempty(Ncirc) 
    Ncirc = 36; %Number of points used in circular inclusion
end 
if nargin < 5 || ~exist('adivd','var') || isempty(adivd) 
    adivd = 1.5; %Period divided by inclusion diameter
end 
if nargin < 6 || ~exist('Nk','var') || isempty(Nk)
    switch kspec                %Number of k-points
        case 'irrfbz'; Nk = 101;  case 'projx';  Nk = [49,11]; 
    end
end 
if nargin < 7 || ~exist('savephi','var') || isempty(savephi)
    savephi = 0; %Whether to save the potential phi (1) or not (0)
end 
if nargin < 9 || ~exist('savename','var') || isempty(savename)
    savename = []; %Unique string identifier may be provided in 'savename'
end 

%Print setup
fprintf('\n|----- SETUP -----|\n')
strf = @(instr,strtype) [repmat(' ',1,14-numel(instr)) instr ' | %' strtype '\n']; %Convenient print string as function
fprintf([strf('lattice type','s'), strf('kspec','s'), strf('Ns','g'), strf('Ncirc','g'), strf('adivd','g'), strf('Nk','s')              , strf('savephi','g')],...
         latticetype             , kspec            , Ns            , Ncirc            , adivd            , strrep(num2str(Nk),'  ','x'), savephi)

%% MESHING & MOMENTA
fprintf('\n|----- MESHING UNIT CELL -----|\n')
theta = linspace(0,2*pi,Ncirc+1); theta = theta(1:end-1); %Inclusion boundary
inclusion = 1/(2*adivd)*[cos(theta);sin(theta)].';        %----||----
cnct = [1:size(inclusion,1);[2:size(inclusion,1),1]].';
[mesh.p,mesh.t] = mesh2d(inclusion,cnct,struct('fun',@(x,y) .025*ones(size(x))),struct('plot',false));
mesh.area = meshArea(mesh.p,mesh.t);
mesh.R = {[0,1],[cosd(60),sind(60)]}
mesh.remesh = mesh;
plotWireMesh(mesh.p,mesh.t)
%mesh = geomPeriodicInclusion(latticetype,inclusion,Ns,@(x,y) .015*250/Ns*ones(size(x)),0);

%Momentum list
switch kspec
    case 'projx'
        [kx,ky] = meshgrid(linspace(0,pi,Nk(1)),linspace(0,1,Nk(2))*2*pi/sqrt(3)); 
        kvec.x = kx.'; kvec.y = ky.'; Nk = Nk(1)*Nk(2); clear kx ky
    case 'irrfbz'
        [kxy,kplot,kmark] = irreducibleFBZ(latticetype,Nk(1),0);
        kvec.x = kxy(:,1); kvec.y = kxy(:,2);  clear kxy
end

%% MODELS 
ef_eV = .48; L = 18e-9*1.5; omegac_eV = @(B) 5.403927588936477e-04*B/ef_eV;
if nargin < 8 || ~exist('models','var') || isempty(models) %Standard, default model
    models(1) = struct('type','isotropic','ef_eV',ef_eV,'L',L,'B',[],'system',[]        ,'approx','intra','keep_eV',[omegac_eV(0),ef_eV]);
    models(2) = struct('type','magneto'  ,'ef_eV',ef_eV,'L',L,'B',4 ,'system','graphene','approx',[]     ,'keep_eV',[omegac_eV(4),ef_eV]); %The lower cutoff is the cyclotron frequency
    models(3) = struct('type','magneto'  ,'ef_eV',ef_eV,'L',L,'B',8 ,'system','graphene','approx',[]     ,'keep_eV',[omegac_eV(8),ef_eV]); %---||---
end
%% ACTUAL CALCULATION AND PLOTTING
bloch.type = '2D'; 

fprintf('|----- LATTICE DISPERSION CALCULATION (%g k-points) -----|\n',Nk)
tick = tic;
for kk2 = size(kvec.x,2):-1:1 %This loop is only relevant in case kspec = 'irrfbz'
    for kk1 = size(kvec.x,1):-1:1
        %--- ACTUAL CALCULATION ---
        bloch.k = [kvec.x(kk1,kk2),kvec.y(kk1,kk2)];
        if savephi == 0 %Only save the eigenenergies
            eneeig_eV_cell = calcEigenEnergiesAny(mesh,models,bloch);
        elseif savephi == 1 %Also save the eigenpotentials
            [eneeig_eV_cell,~,phi_cell] = calcEigenEnergiesAny(mesh,models,bloch);
        end
        for mm = 1:numel(models)
            eneeig_eV{mm}(1:numel(eneeig_eV_cell{mm}),kk1,kk2) = eneeig_eV_cell{mm};
            if savephi == 1
                phieig{mm}{kk1,kk2} = phi_cell{mm};
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
        [rowz,colz] = find(eneeig_eV{mm}(:,:,kk2) == 0);
        eneeig_eV{mm}(sub2ind(size(eneeig_eV{mm}),rowz(rowz~=1),colz(rowz~=1),repmat(kk2,nnz(rowz~=1),1))) = NaN; %This sets any element whose row is larger than 1 and whose value is equal to zero to NaN
    end
end

%% SAVE DATA
clear eneeig_eV_cell phi_cell tick mm kk strf poolobj numProc_cluster numProc_private LASTppn CurMatPool
savedir = '../output/lattice/';
if isempty(savename)
    savename = ['disklattice_' kspec '_Ns' num2str(Ns) '_Ncirc' num2str(Ncirc) '_savephi' num2str(savephi)];
end
savepath = [savedir savename '.mat'];
save(savepath)

%Print final messages
fprintf('\n\n All data saved to   | %s\n\n', savepath)
fprintf('Calculated finished on   | %s\n  ', datestr(now))