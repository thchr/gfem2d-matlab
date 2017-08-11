function runLatticePotentialMap(Ns,Ncirc,Nx)
addpath(genpath('..'))
%% DEFAULTS
if nargin <= 1 || isempty(Ns)
    Ns = 149; %Ns = 230;
      
end
if nargin <= 2 || isempty(Ncirc)
    Ncirc = 64; %Ncirc = 80;
end
disp([Ns,Ncirc])

%% MESH
adivd = 2;
theta = linspace(0,2*pi,Ncirc+1); theta = theta(1:end-1); %Inclusion boundary
inclusion = 1/(2*adivd)*[cos(theta);sin(theta)].';        %----||----
mesh = geomPeriodicInclusion('triangular',inclusion,Ns,@(x,y) .015*250/Ns*ones(size(x)),0);

%% MODEL(S)
ef_eV = .2; L = 400e-9; omegac_eV = @(B) 5.403927588936477e-04*B/ef_eV;
models(1) = struct('type','isotropic','ef_eV',ef_eV,'L',L,'B',0 ,'system',[],        'approx','intra','keep_eV',[0,ef_eV/2]);
models(2) = struct('type','magneto'  ,'ef_eV',ef_eV,'L',L,'B',8 ,'system','graphene','approx',[]     ,'keep_eV',[0,ef_eV/2]); %---||---
%% BLOCH (k-vectors)
%get irrfbz
[fbz.k,fbz.kplot,fbz.kmark] = irreducibleFBZ('triangular',30);
%rotate irrfbz to match the one in our figure
fbz.k = [cosd(60)*fbz.k(:,1) - sind(60)*fbz.k(:,2) , ...
         sind(60)*fbz.k(:,1) + cosd(60)*fbz.k(:,2)];
%pick desired points
klist = { [0,0], fbz.k(fbz.kmark.n(3,1),:) }; 
Ks = exp(1i*[0:5]*2*pi/6)*(klist{2}(:,1) + klist{2}(:,2)*1i); %Get all K points. Lazy way to rotate: (Re,Im) = (x,y)

bloch.type = '2D';
%% RUN

for kk = 1:numel(klist)
    bloch.k = klist{kk};
    [eneeig_eV{kk},rhoeig{kk},phieig{kk},misc{kk}] = calcEigenEnergiesAny(mesh,models,bloch); %Output has form {kk}{mm}(:,modenumber)
end

%% CREATE AN AUXILIARY SQUARE MESH WHICH FILLS IN THE "HOLES" AND EXTENDS
%Create a square plotting mesh
xlims = max(abs(mesh.p(:,1)))*[-1.85,1.85]; 
ylims = max(abs(mesh.p(:,1)))*[-1.65,1.65]; 
if nargin<=3 || isempty(Nx);  %Default
    Nx = 201;
end
if mod(Nx,2) == 0; Nx = Nx + 1; end %Make sure it is an odd integer
Ny = round(diff(ylims)/diff(xlims)*Nx); 
if mod(Ny,2) == 0; Ny = Ny + 1; end %Make sure it is an odd integer

x = linspace(xlims(1),xlims(2),Nx); y = linspace(ylims(1),ylims(2),Ny);
[x,y] = meshgrid(x,y); auxmesh = [x(:),y(:)];

%Find the edge of the unit cell
uc(1,:) = [min(mesh.p(:,1)),min(mesh.p(:,2))];
uc(2,:) = uc(1,:) + mesh.R{1};
uc(3,:) = uc(2,:) + mesh.R{2};
uc(4,:) = uc(3,:) - mesh.R{1};
ucyran = minmax(uc(:,2));


%Create an auxilliary mesh which is folded into the unit cell
inuc = inpolygon(auxmesh(:,1),auxmesh(:,2),uc(:,1),uc(:,2));
for nn = 1:size(auxmesh,1)
    if ~inuc(nn) %If not in unit cell
        breakloop = 1;
        while breakloop
            if auxmesh(nn,2) < ucyran(1)
                auxmesh(nn,:) = auxmesh(nn,:) + mesh.R{2};
            elseif auxmesh(nn,2) > ucyran(2)
                auxmesh(nn,:) = auxmesh(nn,:) - mesh.R{2};
            elseif auxmesh(nn,1) < 0
                auxmesh(nn,1) = auxmesh(nn,1) + mesh.R{1}(:,1);
            elseif auxmesh(nn,1) > 0
                auxmesh(nn,1) = auxmesh(nn,1) - mesh.R{1}(:,1);
            end
            breakloop = ~inpolygon(auxmesh(nn,1),auxmesh(nn,2),uc(:,1),uc(:,2));
        end
    end
end


%% CALCULATE COULOMB OPERATOR ON NEW AUXILIARY MESH
for kk = 1:numel(klist)
    Vaux = calcLatticeCoulombGeneralPoints(mesh,klist{kk},auxmesh);
    for mm = 1:numel(models)
        phiaux{kk}{mm} = Vaux*rhoeig{kk}{mm}; 
    end
end
%% SAVE THE DATA FOR LATER STORAGE
close all; clear cols ans Vaux
save(['../output/lattice/savedPotentialsLatticeTriangular_Ns' num2str(Ns) '_Ncirc' num2str(Ncirc)])

