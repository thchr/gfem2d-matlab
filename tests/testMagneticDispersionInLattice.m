clear all; clc; close all;

%Load necessary subfolders and constants
addpath(genpath('..'));
ConstantsUnits0;

defmesh = 'veryfine';
%Mesh setup
adivd = 2; %Must always be larger than 1! Radius is L/(2*adivd) where a = L
Nk = 25;
latticetype = 'triangular';
meshfun = 'antidot';
plotgraphs = 1;

if exist('defmesh','var')
    switch defmesh
        case 'ultraultrafine' ; Ne = 230;
        case 'ultrafine'; Ne = 149; case 'veryfine' ; Ne = 118;
        case 'fine'     ; Ne = 62;  case 'coarse'   ; Ne = 50;37;
    end
else Ne = 37;
end

if exist('defmesh','var')
    switch defmesh;
        case 'ultraultrafine'; Ncirc = 268;
        case 'ultrafine'; Ncirc = 134; case 'veryfine' ; Ncirc = 66;
        case 'fine'     ; Ncirc = 33;  case 'coarse'   ; Ncirc = 50;16;
    end
else Ncirc = 16;
end

%Inclusion boundary
theta = linspace(0,2*pi,Ncirc+1); theta = theta(1:end-1);
inclusion = 1/(2*adivd)*[cos(theta);sin(theta)].';

%Meshing
blochmesh = geomPeriodicInclusion(latticetype,inclusion,Ne,meshfun,plotgraphs*1);

%Momentum list
[kvec,kplot,kmark] = irreducibleFBZ(latticetype,Nk,plotgraphs*1);

%%
%Material parameters
ef_eV = .2; L = 100e-9;
Blist = [0.01,1:2:100]; eneB_eV = e*Blist*vf^2/(ef_eV*ev2jo)*hbar_eV;

ticEL=tic;
for bb = 1:numel(Blist);
    model = struct('type','magneto','system','graphene','ef_eV',ef_eV,'B',Blist(bb),'L',L); %Model
    [eneeig_eV{bb},eigV{bb}] = calcLatticeEigenEnergies(blochmesh,kvec(kmark.n(3),:),model);
    fprintf('B-FIELD LOOP: %g/%g (%.1f min/%.1f min)\n\n\n',bb,numel(Blist),toc(ticEL)/60,numel(Blist)/bb*toc(ticEL)/60);
end


%%
if exist('plotgraphs','var') && plotgraphs == 1
    
    
    [cols,nams]=flatcolors;
    set_figsize([],14,14); hold on;
    
    plot(eneB_eV/ef_eV,eneB_eV/ef_eV,'--','Color',cols{4})
    for bb = 1:numel(Blist);
        gtr1 = abs(eneeig_eV{bb}/eneB_eV(bb)) >= .999;
        plot(eneB_eV(bb)*ones(size(eneeig_eV{bb}(gtr1))) /ef_eV, real(eneeig_eV{bb}(gtr1))/ef_eV,'.','color',cols{14},'LineWidth',1.5)
        plot(eneB_eV(bb)*ones(size(eneeig_eV{bb}(~gtr1)))/ef_eV, real(eneeig_eV{bb}(~gtr1))/ef_eV,'.','color',cols{8}*.65+cols{4}*.35,'LineWidth',1.5)
    end
    
    hold off
    box on
    set(gca,'Layer','Top','Fontsize',8)
    
    axis tight
    ylim([0,1.05])
    xlim([0,max(eneB_eV)/ef_eV])
   % xlabh = get(gca,'XLabel');
   % set(xlabh,'Position',get(xlabh,'Position') - [0 .1 0])
    xlabel('B-field $\hbar\omega_{\mathrm{c}}/\epsilon_{\mathrm{F}}$','Fontsize',11,'Interpreter','LaTeX')
    ylabel('Frequency $\hbar\omega/\epsilon_{\mathrm{F}}$','Fontsize',11,'Interpreter','LaTeX')
    %   format_ticks(gca,{kmark.symbol{kmark.n(:,2)}},[],[kplot(kmark.n(:,1))],0:1:10,[],[],0.005)
    drawnow;
    
    %Export figure if opened
    export_fig(['../figures/bdisp_hexantidot_diracpoint_' defmesh '_adap'],'-pdf')
    %export_fig(['../figures/magnetodisp_extendedsheet_' defmesh ],'-pdf')
end


%%

Nxy = 25;
x = linspace(-.75,.75,Nxy); y = linspace(-.866/2,.866/2,Nxy); [X,Y] = meshgrid(x,y); r = [X(:),Y(:)];

eW = calcLatticeCoulombOptimized(blochmesh,kvec(kmark.n(3),:));
W = calcLatticeCoulombGeneralPoints(blochmesh,kvec(kmark.n(3),:),r);
R(1,:)=blochmesh.R{1}; R(2,:)=blochmesh.R{2};  

%%

bb=2;

nlist = -1:1; mlist=-1:1;
set_figsize([],50,25);
%eelist = [1,27,35,47:2:63,65,68,69,71];
%eelist = [1,59,75,99:4:131,134,135,138,140];
boti = find(eneeig_eV{bb}/eneB_eV(bb) > 0 & eneeig_eV{bb}/eneB_eV(bb) < 1); topi = find(eneeig_eV{bb}/eneB_eV(bb) > 1); 
eelist = [boti; topi(1:(36-numel(boti)))];
eelist = [boti(end-17:end); topi(1:18)]
for eee=1:36
    ee = eelist(eee);
    
    subaxis(4,9,eee,'M',0.01,'MT',.05,'SH',.01,'SV',.05)
    
    %phi = reshape(W*eigV{bb}(1:size(blochmesh.p,1),ee),size(X));
    phi = blochmesh.remesh.values(eW*eigV{bb}(1:size(blochmesh.p,1),ee));
    hold on
    for n=nlist
        for m=mlist
            Rnm = [n*R(1,1)+m*R(2,1), n*R(1,2)+m*R(2,2)];
            Xnm = blochmesh.remesh.p(:,1)+Rnm(1); Ynm = blochmesh.remesh.p(:,2)+Rnm(2);
            phiexp = phi.*exp(1i*(kvec(kmark.n(3),1)*Xnm + kvec(kmark.n(3),2)*Ynm));
            %imagesc(x+Rnm(1),y+Rnm(2),real(phiexp)); 
            trisurf(blochmesh.remesh.t,Xnm,Ynm,real(phiexp),'EdgeColor','none');
            plot(inclusion(:,1)+Rnm(1),inclusion(:,2)+Rnm(2),'--','color',cols{3});
        end
    end
    plot([-.75,.25,.75,-.25,-.75],[-1,-1,1,1,-1]*R(2,2)/2,'--','color',cols{3}); 
    hold off
    set(gca,'YDir','normal')
    
    cvals=sort(abs(real(phiexp(:))),'descend');
    %caxis([-1,1]*cvals(20))
    shading interp;
    colormap('bluewhitered_mod'); 
    freezeColors;
    view(2);
    axis equal off;
    xlim([-1,1]*1.35); ylim([-1,1]*1.2)
    title(['\omega/\omega_c = ' num2str(real(eneeig_eV{bb}(ee)/eneB_eV(bb)),2)]);
    drawnow; %pause(.25)
end