clear all; clc; close all;

addpath(genpath('..'))

ConstantsUnits0;
N = 50;
geom = 'disk';
switch geom
    case 'disk'
        [~,p,t] = evalc('geomDisk(N,7.5/N,1)');
    case 'triangle'
        [~,p,t] = evalc('geomEquiTriangle(1,.0625,1)');
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

ef_eV = .15;
gam_eV = 24e-3;
L = 100e-9;

%%
Blist = [0.01,1:1:40]; eneB_eV = e*Blist*vf^2/(ef_eV*ev2jo)*hbar_eV;

ticEL=tic;
for bb = 1:numel(Blist);
    model = struct('type','magneto','system','graphene','ef_eV',ef_eV,'B',Blist(bb),'L',L,'lossc',gam_eV/eneB_eV(bb)); %Model
    [eneeig_eV{bb},eigV{bb}] = calcEigenEnergies(struct('p',p,'t',t),model);
    fprintf('B-FIELD LOOP: %g/%g (%.1f min/%.1f min)\n\n\n',bb,numel(Blist),toc(ticEL)/60,numel(Blist)/bb*toc(ticEL)/60);
end

model0 = struct('type','isotropic','approx','intra','ef_eV',ef_eV,'L',L);
[eneeig0_eV,eigV0] = calcEigenEnergies(struct('p',p,'t',t),model0);

%%
    [cols,nams]=flatcolors;
    set_figsize([],14,14); hold on;
    
    plot(eneB_eV/ef_eV,eneB_eV/ef_eV,'--','Color',cols{4})
    for bb = 1:numel(Blist);
        gtr1 = abs(real(eneeig_eV{bb})/eneB_eV(bb)) >= .999; %| abs(imag(eneeig_eV{bb}/eneB_eV(bb)))<.5*gam_eV/eneB_eV(bb);
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
    drawnow
    break
%%

[~,V] = evalc('calcCoulomb(p,t)');
%%
bb=22;
ylims = [0,2500000];
eelist=[];
for ee=1:numel(eneeig_eV{bb})
    if real(eneeig_eV{bb}(ee))/eneB_eV(bb) > ylims(1) && real(eneeig_eV{bb}(ee))/eneB_eV(bb) < ylims(2)
        real(eneeig_eV{bb}(ee))
        eelist = [eelist;ee];
    end
end
disp(eneB_eV(bb)/ef_eV)
eelists = eelist([1,5,10,13,16:end])
eelists = eelist([1,5,10,12,13,18:end])
try; close(4); end
set_figsize(4,42,25);
for eee=1:36
    ee = eelists(eee);
    subaxis(4,7,eee,'M',0.01,'MT',.05,'SH',.01,'SV',.05)
    trisurf(t,p(:,1),p(:,2),real(V*eigV{bb}(:,ee)),'EdgeColor','none'); hold on
    shading interp;
    switch geom; 
        case 'ring'
            plot(nodesout(:,1),nodesout(:,2),'--','Color',cols{3})
            plot(nodesin(:,1),nodesin(:,2),'--','Color',cols{3}); hold off
        case 'disk'
            plot(cos(linspace(0,2*pi,100)),sin(linspace(0,2*pi,100)),'--','Color',cols{3}); hold off
    end
    colormap('bluewhitered_mod'); freezeColors;
    view(2);
    axis equal off;
    title(['\omega/\omega_c = ' num2str(real(eneeig_eV{bb}(ee)/eneB_eV(bb)),3)],'Fontsize',10); drawnow; 
end

export_fig(['../figures/potentials_eig_magneto ' geom '_enecdivef_' strrep(num2str(eneB_eV(bb)/ef_eV,2),'.','p')],'-jpeg')

%% gradient map
bb=22;
rc = meshCentroid(p,t)
try; close(5); end
set_figsize(5,42,25);
eelists = eelist([1,5,10,13,16:end])
eelists = eelist([1,5,10,12,13,18:end])
for eee= [1,1:28];
    ee = eelists(eee);
    E = -calcGradient(p,t,V*eigV{bb}(:,ee)); 
    sigma = conducMagnetoClassical(eneeig_eV{bb}(ee),ef_eV,Blist(bb),0);
    J = (sigma*E.').';
    
    
    subaxis(4,7,eee,'M',0.01,'MT',.05,'SH',.01,'SV',.05)
    quiverc(rc(:,1),rc(:,2),imag(J(:,1)),imag(J(:,2))); hold on
    switch geom; 
        case 'ring'
            plot(nodesout(:,1),nodesout(:,2),'--','Color',cols{3})
            plot(nodesin(:,1),nodesin(:,2),'--','Color',cols{3}); hold off
        case 'disk'
            plot(cos(linspace(0,2*pi,100)),sin(linspace(0,2*pi,100)),'--','Color',cols{3}); hold off
    end
    colormap(bluewhitered_mod); freezeColors;
    view(2);
    axis equal off;
   % title([num2str(ee) ', ' num2str(real(eneeig_eV{bb}(ee)/eneB_eV(bb)),3)]); 
    title(['\omega/\omega_c = ' num2str(real(eneeig_eV{bb}(ee)/eneB_eV(bb)),3)],'Fontsize',10,'Color','w');
    drawnow; 
    set(gcf,'Color',[0,0,0])
end

export_fig(['../figures/current_eig_magneto ' geom '_enecdivef_' strrep(num2str(eneB_eV(bb)/ef_eV,2),'.','p')],'-pdf')