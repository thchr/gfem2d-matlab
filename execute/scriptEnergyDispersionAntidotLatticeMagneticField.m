clear all; clc; close all;

%Load necessary subfolders and constants
addpath(genpath('..'));
ConstantsUnits0; 

defmesh = 'fine';
%Check input
if ~exist('adivd','var') || isempty(adivd)
	adivd = 1.5; %Must always be larger than 1! Radius is L/(2*adivd) where a = L
end
if ~exist('Ne','var') || isempty(Ne)
    if exist('defmesh','var')
        switch defmesh
            case 'ultraultrafine' ; Ne = 230; 
            case 'ultrafine'; Ne = 149; case 'veryfine' ; Ne = 118;
            case 'fine'     ; Ne = 62;  case 'coarse'   ; Ne = 37;
        end
    else Ne = 37;
    end
end
if ~exist('Ncirc','var') || isempty(Ncirc)
    if exist('defmesh','var')
        switch defmesh; 
            case 'ultraultrafine'; Ncirc = 268;
            case 'ultrafine'; Ncirc = 134; case 'veryfine' ; Ncirc = 66;
            case 'fine'     ; Ncirc = 33;  case 'coarse'   ; Ncirc = 16;
        end
    else Ncirc = 16;
    end
end
if ~exist('Nk','var') || isempty(Nk)
	Nk = 35;
end
if ~exist('latticetype','var') || isempty(latticetype)
	latticetype = 'triangular';
end
if ~exist('meshfun','var') || isempty(meshfun)
	meshfun = 'const';
end
if ~exist('plotgraphs','var') || isempty(plotgraphs)
    plotgraphs = 1;
end



%Inclusion boundary
theta = linspace(0,2*pi,Ncirc+1); theta = theta(1:end-1);
inclusion = 1/(2*adivd)*[cos(theta);sin(theta)].';

%Meshing
blochmesh = geomPeriodicInclusion(latticetype,inclusion,Ne,meshfun,plotgraphs*1);

%Momentum list
[kvec,kplot,kmark] = irreducibleFBZ(latticetype,Nk,plotgraphs*1);

%Model 
B = 10; ef_eV = .2; L = 100e-9; eneB_eV = e*B*vf^2/(ef_eV*ev2jo)*hbar_eV;
model = struct('type','magneto','system','graphene','ef_eV',ef_eV,'B',B,'L',L);
%model = struct('type','isotropic','approx','intra','ef_eV',ef_eV,'L',L,'keepzerosols',1);


%%
%parpool;
ticEL=tic;
for kk = 1:size(kvec,1);
    [eneeig_eV{kk},eigV{kk}] = calcLatticeEigenEnergies(blochmesh,kvec(kk,:),model);
    fprintf('DISPERSION LOOP: %g/%g (%.1f min/%.1f min)\n\n\n',kk,size(kvec,1),toc(ticEL)/60,size(kvec,1)/kk*toc(ticEL)/60);
end


%%
if exist('plotgraphs','var') && plotgraphs == 1

	
	[cols,nams]=flatcolors;
	set_figsize([],14,14); hold on;

	for mm = 2:size(kmark.n,1)-1
		plot(kplot(kmark.n(mm,1))*[1,1],[0,1]*20,'-','Color',[1,1,1]*.25,'LineWidth',.45)
	end
	for kk = 1:size(kvec,1)
        gtr1 = abs(eneeig_eV{kk}/eneB_eV) >= .999;
        plot(kplot(kk)*ones(size(eneeig_eV{kk}(gtr1))), real(eneeig_eV{kk}(gtr1))/eneB_eV,'.','color',cols{2},'LineWidth',1.5)
        plot(kplot(kk)*ones(size(eneeig_eV{kk}(~gtr1))), real(eneeig_eV{kk}(~gtr1))/eneB_eV,'.','color',cols{2},'LineWidth',1.5)
    end
	hold off
	box on
	set(gca,'Layer','Top','Fontsize',8)
	
	axis tight
	ylim([0,1]*ef_eV/eneB_eV)
	xlim([0,max(kplot)])
	xlabh = get(gca,'XLabel');
	set(xlabh,'Position',get(xlabh,'Position') - [0 .1 0])
	xlabel('Momentum','Fontsize',10)
	ylabel('Frequency \omega/\omega_c','Fontsize',10)
    format_ticks(gca,{kmark.symbol{kmark.n(:,2)}},[],[kplot(kmark.n(:,1))],0:1:10,[],[],0.005)
    drawnow;
    
    %Export figure if opened
	%export_fig(['../figures/magnetodisp_hexantidot_' defmesh '_adapmesh' '_spurious' ],'-pdf')
    %export_fig(['../figures/magnetodisp_extendedsheet_' defmesh ],'-pdf')
end

%%
break

kk = 1;kmark.n(3,1);

Nxy = 50;
x = linspace(-.75,.75,Nxy); y = linspace(-.866/2,.866/2,Nxy); [X,Y] = meshgrid(x,y); r = [X(:),Y(:)];

%W = calcLatticeCoulombOptimized(blochmesh,kvec(kk,:));
W = calcLatticeCoulombGeneralPoints(blochmesh,kvec(kk,:),r);
R(1,:)=blochmesh.R{1}; R(2,:)=blochmesh.R{2};  

%%


nlist = -2:2; mlist=-2:2;
set_figsize([],50,25);
%eelist = [1,27,35,47:2:63,65,68,69,71];
%eelist = [1,59,75,99:4:131,134,135,138,140];
eelist = [4,26,50,84,90,96,102,108,113,119,125,131,134,135,138,140];
boti = find(eneeig_eV{kk}/eneB_eV > 0 & eneeig_eV{kk}/eneB_eV < 1); topi = find(eneeig_eV{kk}/eneB_eV > 1); 
eelist = [boti; topi(1:(36-numel(boti)))];
for eee=1:36
    ee = eelist(eee);
    
    subaxis(4,9,eee,'M',0.01,'MT',.05,'SH',.01,'SV',.05)
    %trisurf(blochmesh.remesh.t,blochmesh.remesh.p(:,1),blochmesh.remesh.p(:,2),real(blochmesh.remesh.values(((eigV{kk}(1:size(blochmesh.p,1),ee))))),'EdgeColor','none');
    phi = reshape(W*eigV{kk}(1:size(blochmesh.p,1),ee),size(X));
    hold on
    for n=nlist
        for m=mlist
            Rnm = [n*R(1,1)+m*R(2,1), n*R(1,2)+m*R(2,2)];
            Xnm = X+Rnm(1); Ynm = Y+Rnm(2);
            phiexp = phi.*exp(1i*(kvec(kk,1)*Xnm + kvec(kk,2)*Ynm));
            imagesc(x+Rnm(1),y+Rnm(2),real(phiexp)); 
            plot(inclusion(:,1)+Rnm(1),inclusion(:,2)+Rnm(2),'--','color',cols{3});
        end
    end
    plot([-.75,.25,.75,-.25,-.75],[-1,-1,1,1,-1]*R(2,2)/2,'--','color',cols{3}); 
    hold off
    set(gca,'YDir','normal')
    
    cvals=sort(abs(real(phiexp(:))),'descend');
    caxis([-1,1]*cvals(20))
    %shading interp;
    colormap('bluewhitered_mod'); 
    freezeColors;
    view(2);
    axis equal off;
    xlim([-1,1]*1.35); ylim([-1,1]*1.2)
    title(['\omega/\omega_c = ' num2str(real(eneeig_eV{kk}(ee)/eneB_eV),2)]);
    drawnow;% pause(.25)
end
%export_fig(['../figures/surfs_selecteigenpotentials_hexantidot_' defmesh '_adapmesh'],'-pdf')
