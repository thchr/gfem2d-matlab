clc; clear all; close all;
cols = flatcolors; 
%% MESH
Ns = 230; Ncirc = 80; 
%Ns = 149; Ncirc = 64; 
adivd = 2;
theta = linspace(0,2*pi,Ncirc+1); theta = theta(1:end-1); %Inclusion boundary
inclusion = 1/(2*adivd)*[cos(theta);sin(theta)].';        %----||----
mesh = geomPeriodicInclusion('triangular',inclusion,Ns,@(x,y) .015*250/Ns*ones(size(x)),1);

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
%Illustrate selected k-points in BZ
set_figsize(2,10,10)
plot(real(Ks([1:end,1])),imag(Ks([1:end,1])),':k'); hold on
plot(fbz.k(:,1),fbz.k(:,2),'-','color',cols{14},'LineWidth',2.5);
for kk = 1:numel(klist); 
    plot(klist{kk}(:,1),klist{kk}(:,2),'o','color',cols{1},'LineWidth',2,'markerfacecolor','w'); 
end
axis equal; hold off;
bloch.type = '2D';
%% RUN

for kk = 1:numel(klist)
    bloch.k = klist{kk};
    [eneeig_eV{kk},rhoeig{kk},phieig{kk},misc{kk}] = calcEigenEnergiesAny(mesh,models,bloch); %Output has form {kk}{mm}(:,modenumber)
end

%% CREATE AN AUXILIARY SQUARE MESH WHICH FILLS IN THE "HOLES" AND EXTENDS
%Create a square plotting mesh
Nsq = 250;
xlims = max(abs(mesh.p(:,1)))*[-1.85,1.85]; 
ylims = max(abs(mesh.p(:,1)))*[-1.65,1.65]; 
x = linspace(xlims(1),xlims(2),Nsq); y = linspace(ylims(1),ylims(2),Nsq);
[x,y] = meshgrid(x,y); auxmesh = [x(:),y(:)];

%Find the edge of the unit cell
uc(1,:) = [min(mesh.p(:,1)),min(mesh.p(:,2))];
uc(2,:) = uc(1,:) + mesh.R{1};
uc(3,:) = uc(2,:) + mesh.R{2};
uc(4,:) = uc(3,:) - mesh.R{1};
ucyran = minmax(uc(:,2));
set_figsize(3,10,10)
plot(uc([1:end,1],1),uc([1:end,1],2),'-','color',cols{4})
hold on; plot(mesh.p(:,1),mesh.p(:,2),'+','color',cols{1}); 

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
plot(auxmesh(:,1),auxmesh(:,2),'.','color',cols{14}); axis equal; hold off;

%% CALCULATE COULOMB OPERATOR ON NEW AUXILIARY MESH
for kk = 1:numel(klist)
    Vaux = calcLatticeCoulombGeneralPoints(mesh,klist{kk},auxmesh);
    for mm = 1:numel(models)
        phiaux{kk}{mm} = Vaux*rhoeig{kk}{mm}; 
    end
end
%% SAVE THE DATA FOR LATER STORAGE
close all; clear cols ans
save(['savedPotentialsLatticeTriangular_Ns' num2str(Ns) '_Ncirc' num2str(Ncirc)])

%% PLOT THE POTENTIALS ON THE GRIDDED MESH
%load('savedPotentialsTriangular_Ns149_Ncirc64')
mmchoice = 2; 
kkchoice = 2; 
nnchoice = 2; %194 (k=0), [216, 215, (k=K)]
plotvals = real(reshape(phiaux{kkchoice}{mmchoice}(:,nnchoice),size(x)));
set_figsize(4,10,10)
%draw the potential
contourf(x,y,plotvals,25,'EdgeColor','none');
%surf(x,y,plotvals); shading interp; 
colormap(bluewhitered_mod(512)); 
%colormap(flipud(morgenstemning)); 
caxis(max(abs(plotvals(:)))*[-1,1])
hold on 

%draw the unit cell
plot(uc([1:end,1],1),uc([1:end,1],2),'-','color',cols{4})

%draw boundary (i.e. circles)
circstyle = {'Curvature',[1,1],'edgecolor',cols{8}*.25+cols{4}*.75,'linestyle',':'};
rectangle('Position',[-.5,-.5,1,1]*adivd/4,circstyle{:}); %One can draw circles with the rectangle command
for cc = 0:5
    rectangle('Position',[cosd(60*cc),sind(60*cc),0,0]*1+[-.5,-.5,1,1]*adivd/4,circstyle{:});
end
hold off
%axis off; 
axis equal tight
xlim(xlims); ylim(ylims)

%% IDENTIFY THE k=0 FINITE B-FIELD MODE
enecase = eneeig_eV{kkchoice}{mmchoice}; 
ylimsearch = [1e-5,10*omegac_eV(models(2).B)];

set_figsize(5,25,25)
plot(1:numel(enecase),enecase,'.')
ylim(ylimsearch); 
xlim([min(find(enecase>ylimsearch(1))),max(find(enecase<ylimsearch(2)))])
%use this plot to find the non-monotonically increasing mode