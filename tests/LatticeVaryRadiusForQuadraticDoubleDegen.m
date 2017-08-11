clear; clc; close all;
Ns=[]; Ncirc = []; Nk = []; savephi = 0;
adivdL = 2.5:-.1:1.5;
ef_eV = .4; L = 400e-9;
models = struct('type','isotropic','ef_eV',ef_eV,'L',L,...
                'approx','intra','keep_eV',[0,1.5*ef_eV]);
for adivd = adivdL
    addstr = ['_adivd' num2str(adivd,'%.2f')];
    runLatticeMagneto('triangular','irrfbz',Ns,Ncirc,adivd,Nk,savephi,models,addstr)
end

%%
MKi = 37:58;    %path from M->K
KGi = 58:101;   %path from K->G
GMi = 1:37;     %path from G->M
shifti = [MKi,KGi(2:end),GMi(2:end)];  %path from M->K->G->M
shifti = fliplr(shifti);               %path from M->G->K->M
newtickmarks = {'\bfM', '\bf\Gamma', '\bfK', '\bfM'};
%%
theta0 = linspace(0,2*pi,101); theta0 = theta0(1:end-1);
xc0 = cos(theta0); yc0 = sin(theta0); 

cols = flatcolors(); ConstantsUnits0;
close all;
for aa = 1:numel(adivdL)
    loadstr = ['../output/lattice/lattice_irrfbz_Ns70_Ncirc36_savephi0__adivd' num2str(adivdL(aa),'%.2f') '.mat'];
    load(loadstr)
    
    kvec.x = kvec.x(shifti,:);
    kvec.y = kvec.y(shifti,:);
    newticks = kplot([1,GMi(end),numel([KGi,GMi(2:end)]),numel(shifti)]);
    
    set_figsize(aa,5.75,10);

    subaxis(1,1,1,'MR',0.04,'ML',.2,'MT',.015)
    
    for mm = 2:size(kmark.n,1)-1
        plot(newticks(mm)*[1,1],[0,1]*120,'-','Color',cols{8}*.5+cols{4}*.5,'LineWidth',.4); hold on
    end
    plot(kplot,real(eneeig_eV{1}(:,shifti))*eV2THz,'-','color',cols{14},'LineWidth',1.25)
    
    
    axis tight
    xlim([0,max(kplot)])
    ylim([0,ef_eV*.413567]*eV2THz)

    set(gca,'Fontsize',7,'Layer','Top','LineWidth',.4,'TickLength',[1,1]*0.015,'TickDir','out','XTick',newticks,'XTicklabel',newtickmarks)
    ylabel('Frequency (THz)','Fontsize',9)
    xlabel('Wave vector \bfk','Fontsize',9)
    
    
    %Fine tuning by hand is necessary here
    handaxes2 = axes('Position', [0.575 0.065 0.35 0.35]);
    %Plot the mesh interior
    rectangle('Position',[-1.51,-1.6,3.02,2.75],'facecolor','w','edgecolor','none')
    patch([-1.5,.5,1.5,-.5],[-1,-1,1,1]*sqrt(3)/2,cols{19}*.9)
    patch(xc0/adivd,yc0/adivd,'w');
    set(handaxes2, 'box', 'off','color','none'); axis off
    text(-1.3,-1.25,['\itd\rm/\ita\rm = 1/' num2str((adivd),'%.1f')],'Fontsize',8);
    drawnow
    axis equal
end