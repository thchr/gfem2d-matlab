clear all; close all; clc; 
ConstantsUnits0; 

%% SET UP A NUMBER OF MODELS
ef_eV = .2; L = 400e-9;
models(1) =     struct('type','isotropic','ef_eV',ef_eV,'L',L,'F',0 ,'approx','intra','keep_eV',[0,ef_eV]);
models(end+1) = struct('type','berry'    ,'ef_eV',ef_eV,'L',L,'F',.1,'approx',[]     ,'keep_eV',[0,ef_eV]); 
models(end+1) = struct('type','berry'    ,'ef_eV',ef_eV,'L',L,'F',.2,'approx',[]     ,'keep_eV',[0,ef_eV]); 
models(end+1) = struct('type','berry'    ,'ef_eV',ef_eV,'L',L,'F',.3,'approx',[]     ,'keep_eV',[0,ef_eV]);
models(end+1) = struct('type','berry'    ,'ef_eV',ef_eV,'L',L,'F',.4,'approx',[]     ,'keep_eV',[0,ef_eV]);
models(end+1) = struct('type','berry'    ,'ef_eV',ef_eV,'L',L,'F',.5,'approx',[]     ,'keep_eV',[0,ef_eV]);
models(end+1) = struct('type','isotropic','ef_eV',ef_eV*2,'L',L/2,'F',0 ,'approx','intra','keep_eV',[0,ef_eV*2]);
models(end+1) = struct('type','berry'    ,'ef_eV',ef_eV*2,'L',L/2,'F',.1,'approx',[]     ,'keep_eV',[0,ef_eV*2]); 
models(end+1) = struct('type','berry'    ,'ef_eV',ef_eV*2,'L',L/2,'F',.2,'approx',[]     ,'keep_eV',[0,ef_eV*2]); 
models(end+1) = struct('type','berry'    ,'ef_eV',ef_eV*2,'L',L/2,'F',.3,'approx',[]     ,'keep_eV',[0,ef_eV*2]);
models(end+1) = struct('type','berry'    ,'ef_eV',ef_eV*2,'L',L/2,'F',.4,'approx',[]     ,'keep_eV',[0,ef_eV*2]);
models(end+1) = struct('type','berry'    ,'ef_eV',ef_eV*2,'L',L/2,'F',.5,'approx',[]     ,'keep_eV',[0,ef_eV*2]);

%% EXECUTE BAND DIAGRAM CALCULATION FOR REQUESTED MODELS
%Ns = 118; Ncirc = 50; 
Ns = []; Ncirc = []; 
Nk = 101;35; 
savepath = runLatticeMagneto('triangular','irrfbz',Ns,Ncirc,[],Nk,[],models,'berry');

%% LOAD DATA JUST OBTAINED (saved to harddisk, not in memory)
load(savepath);

%% PLOT DATA 
cols = flatcolors; 
%indices to enter into bands{mm}.ene_eV{bb}(...)
MKi = kmark.n(2,1):kmark.n(3,1);    %path from M->K
KGi = kmark.n(3,1):kmark.n(4,1);   %path from K->G
GMi = kmark.n(1,1):kmark.n(2,1);     %path from G->M
shifti = [MKi,KGi(2:end),GMi(2:end)];  %path from M->K->G->M
shifti = fliplr(shifti);

xticks = kplot([1,GMi(end),numel([KGi,GMi(2:end)]),numel(shifti)]); 


%plotting
try; close(2); end
set_figsize(2,3.5*numel(models)/2+3,9*1.85)
for mm = 1:numel(models)
    ylims = [0,0.5*models(mm).ef_eV]*eV2THz; 
    yticks = 0.516958439800970*models(mm).ef_eV*eV2THz*(0:4)/5;
    
    subaxis(2,numel(models)/2,mm,'SH',.01,'MR',.01,'ML',.05,'SV',0.014,'MT',.02,'MB',.06)
    %plotting data
    for nn = 2:numel(xticks)-1
    plot(xticks(nn)*[1,1],ylims,'-','color',cols{8},'linewidth',.4); hold on
    end
    plot(kplot,real(eneeig_eV{mm}(:,shifti))*eV2THz,'color',cols{14},'linewidth',.8);
    hold off
    

	%annotations, clean up, etc
    if ~any(mm == [1,numel(models)/2+1])
        set(gca,'YTickLabel',{''})
        berrytxt = []; 
    else
        ylabel('Frequency (THz)','Fontsize',9)
        berrytxt = 'Net Berry flux, ';
    end
    
    text(kplot(50),ylims(2)*.075,{['\itE\rm_F = ' num2str(models(mm).ef_eV) ' eV'],...
                                  ['\ita\rm = ' num2str(models(mm).L*1e9) ' nm']},...     
         'BackgroundColor','w','Fontsize',9)
	text(kplot(end),ylims(2),[' ' berrytxt '\itF\rm = ' num2str(models(mm).F) ' '],...
         'BackgroundColor','k','Color','w','Fontsize',9,'VerticalAlignment','top','HorizontalAlignment','right','Margin',0.1)
     
    ylim(ylims)
    xlim(minmax(kplot))
    set(gca,'XTick',xticks,'XTickLabel',{'M','G','K','M'},'Fontsize',7,...
        'YTick',yticks,'LineWidth',.4,'TickLength',[.02,.02])
    if mm >= numel(models)/2+1
        xlabel('Wave vector \bfk','Fontsize',9)
    else
        set(gca,'XTickLabel','')
    end
end

%save figure
%export_fig('../figures/lattice/berryflux_triangular','-pdf')
