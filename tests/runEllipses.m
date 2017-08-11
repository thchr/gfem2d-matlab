clear all; close all; clc;

addpath(genpath('..'))

ConstantsUnits0;
N = 100;

alist = linspace(sqrt(1/pi),2,7);
ellipse = cell(numel(alist));
for aa = 1:numel(alist);
    a = alist(aa); b = 1/(pi*a); blist(aa) = b;
    
    theta = linspace(0,2*pi,N+1).'; theta = theta(1:end-1) + (theta(2)-theta(1))/2;
    nodes = [a*cos(theta),b*sin(theta)];
    [~,mesh.p,mesh.t] = evalc('geomPolygon(nodes,5/N,1)');
    
    [alpha,eigstruct] = calcPolarizability(mesh,[],1);
    
    %Store the data in a cell array of structures
    ellipse{aa}.a = a; ellipse{aa}.b = b;
    ellipse{aa}.bndry = [nodes; nodes(1,:)];
    ellipse{aa}.mesh = mesh;
    ellipse{aa}.alpha = alpha;
    ellipse{aa}.eigstruct = eigstruct;
    ellipse{aa}.N = N;
end
%%
ef_eV = .4;
gam_eV = 12e-3;
L = 20e-9;
ene_eV = linspace(0,2,2500)*ef_eV; omega = ene_eV/hbar_eV;
sigma =  conducLRA(ene_eV,ef_eV,gam_eV,300*kb_eV);
zeta = 2i*eps0*omega*L./sigma;

%%
cols=flatcolors; colx = cols{1}; coly = cols{14}*.75+cols{24}*.25;
set_figsize([],17,25);
for aa = 1:numel(alist)
    subaxis(numel(alist),1,aa,'SV',0.0025,'MT',0.025,'MR',.025,'MB',.05)
    Qabsx = omega/c.*imag(ellipse{aa}.alpha.x(zeta,L))/L^2; 
    Qabsy = omega/c.*imag(ellipse{aa}.alpha.y(zeta,L))/L^2; 

    [h1,h2]=ShadedDualArea(ene_eV,Qabsx,Qabsy,0,{colx coly},{0.85,.9}); hold on
    
    if aa ~= numel(alist)
        set(gca,'XTickLabel',{''})
    else
        xlabel('Energy [eV]','Fontsize',10)
    end
    if aa == 1
        legend([h1,h2],{'\itx\rm-polarized','\ity\rm-polarized'},'Fontsize',9,'Location','NorthWest'); legend boxoff;
    end
    if aa == 4
        ylabel('Absorption efficiency  \sigma_{abs}/area','Fontsize',10)
    end
    set(gca,'Fontsize',7,'LineWidth',.15)
    ylim([0,1]); set(gca,'YTick',0:.3:.9)

    %Inset of mesh geometry (fine tuning by hand is necessary for position)
    handaxes2 = axes('Position', [0.845 .875-.132*(aa-1) 0.115 0.15]);
    %Plot the polygon boundary
    htp=patch(ellipse{aa}.bndry(:,1),ellipse{aa}.bndry(:,2),cols{19}*.9,'EdgeColor','none');
    text(0,-ellipse{aa}.b-.85/sqrt(pi),['\ita\rm = ' num2str(ellipse{aa}.a/ellipse{aa}.b,2) '\itb\rm'],'Fontsize',8,'HorizontalAlignment','c')
    axis equal off
    xlim([-1.025,1.025]*max(alist))
    ylim([-1.025,1.025]*1/sqrt(pi))
    set(handaxes2, 'box', 'off','color','none');
    drawnow
    
end

hold off

%Export figure 
%export_fig(['AbsorptionEfficiency_Ellipse_EqualArea400nm2_Ef400meV_gam12meV_N' num2str(N)],'-pdf')


%%
AA = 2;
plotVertexData(ellipse{AA}.mesh.p,ellipse{AA}.mesh.t,ellipse{AA}.eigstruct.phi(:,7))