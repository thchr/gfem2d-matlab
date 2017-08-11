clear all; close all; clc;

addpath(genpath('..'))

ConstantsUnits0;
N = 75;
geom = 'disk';
switch geom
    case 'disk'
        [~,p,t] = evalc('geomDisk(N,7.5/N,1)');
    case 'ellipse'
        a = 1.1; b = .9;
        theta = linspace(0,2*pi,N+1).'; theta = theta(1:end-1) + (theta(2)-theta(1))/2;
        nodes = [a*cos(theta),b*sin(theta)];
        [~,p,t] = evalc('geomPolygon(nodes,5/N,1)');
    case 'ring'
        rin = .15;
        thetaout = linspace(0,2*pi,N+1).'; thetaout = thetaout(1:end-1) + (thetaout(2)-thetaout(1))/2;
        thetain = linspace(0,2*pi,round(N*sqrt(rin))+1).'; thetain = thetain(1:end-1) + (thetain(2)-thetain(1))/2;
        nodesout = [cos(thetaout),sin(thetaout)]; nodesin = rin*[cos(thetain),sin(thetain)];
        nodes = {nodesout ,nodesin};
        [~,p,t] = evalc('geomPolygon(nodes,5/N,1)');
end


%
ef_eV = .2;
gam_eV = 6e-3;
L = 100e-9;
ene_eV = linspace(.15,0.55,50)*ef_eV; omega = ene_eV/hbar_eV;
%sigma =  conducLRA(ene_eV,ef_eV,gam_eV,300*kb_eV);
B = 10;
eneB_eV = e*B*vf^2/(ef_eV*ev2jo)*hbar_eV;


sigmalra = @(omega) conducIntra(omega*hbar_eV,ef_eV,gam_eV);
sigma = @(omega) conducMagnetoClassical(omega*hbar_eV,ef_eV,B,gam_eV);
%sigma = @(omega) calcMagnetoClassical(omega*hbar_eV,ef_eV,B,gam_eV);
zeta = 2i*eps0*omega*L./sigmalra(omega);
f = @(omega) sigma(omega)./sigmalra(omega);

%%


[Q_abs] = calcAbsorption(p,t,omega,L,0,sigmalra,f);
[alpha0,eigstruct] = calcPolarizability(struct('p',p,'t',t));%round(size(mesh.p,1)/2));
Q_abs0 = omega./c.*imag(alpha0.x(zeta,L))/(pi*L^2);
%%

model = struct('type','magneto','system','graphene','ef_eV',ef_eV,'B',B,'L',L);
[eneeig_eV,eigV] = calcEigenEnergies(struct('p',p,'t',t),model);

model0 = struct('type','isotropic','approx','intra','ef_eV',ef_eV,'L',L);
[eneeig0_eV,eigV0] = calcEigenEnergies(struct('p',p,'t',t),model0);

%%
try; close(2); end
cols=flatcolors;
set_figsize(2,25,15);

plot(ene_eV,Q_abs,'-','color',cols{14}); hold on;
plot(ene_eV,Q_abs0,'--','color',cols{1})

for n=1:numel(eneeig_eV)
    if abs(eneeig_eV(n)) > min(ene_eV) && abs(eneeig_eV(n)) < max(ene_eV)
        plot(eneeig_eV(n)*[1,1],[0,max(Q_abs)],'-','color',cols{14})
        text(eneeig_eV(n),max(Q_abs)*1.05,num2str(n),'Fontsize',8,'Color',cols{14},'HorizontalAlignment','Center','VerticalAlignment','Bottom')
    end
end

for n = 1:numel(eneeig0_eV)
    plot(eneeig0_eV(n)*[1,1],[0,max(Q_abs0)],'--','color',cols{1})
end

hold off
xlim(minmax(ene_eV))
%%
ylims = [0,7];
try; close(3); end
eelist=[];
set_figsize(3,25,15); hold on
for ee=1:numel(eneeig_eV)
    if real(eneeig_eV(ee))/eneB_eV > ylims(1) && real(eneeig_eV(ee))/eneB_eV < ylims(2)
        eelist = [eelist;ee];
        plot([0,1],[1,1]*real(eneeig_eV(ee))/eneB_eV,'-','color',cols{2})
    end
end
for ee=1:numel(eneeig0_eV)
if real(eneeig0_eV(ee))/eneB_eV > ylims(1) && real(eneeig0_eV(ee))/eneB_eV < ylims(2)
        plot([0,1],[1,1]*real(eneeig0_eV(ee))/eneB_eV,'--','color',cols{1})
    end
end
%ylim(ylims)
%%
%Sort according (small to large)
%[lambdat2,I]=sort(real(lambdat)); eigVB2 = eigVB(:,I);

%plotVertexData(p,t,real(eigVB(:,447)))
%plotVertexData(p,t,real(eigVB(:,448)))
%break

[~,V] = evalc('calcCoulomb(p,t)');
%%
try; close(4); end
set_figsize(4,50,25);
for eee=1:36
    ee = eelist(end-10-eee);
    subaxis(4,9,eee,'M',0.01,'MT',.05,'SH',.01,'SV',.05)
    trisurf(t,p(:,1),p(:,2),real(V*eigV(:,ee)),'EdgeColor','none'); hold on
    shading interp;
    switch geom; 
        case 'ring'
            plot(nodesout(:,1),nodesout(:,2),'--','Color',cols{3})
            plot(nodesin(:,1),nodesin(:,2),'--','Color',cols{3}); hold off
        case 'disk'
     %       plot(nodes(:,1),nodes(:,2),'--','Color',cols{3}); hold off
    end
    colormap('bluewhitered_mod'); freezeColors;
    view(2);
    axis equal off;
    title([num2str(ee) ', ' num2str(real(eneeig_eV(ee)/eneB_eV),3)]); drawnow; pause(.1)
end