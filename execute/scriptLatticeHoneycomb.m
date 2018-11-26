clear all; close all; clc;
ConstantsUnits0;

%% EXECUTE BAND DIAGRAM CALCULATION FOR REQUESTED MODELS
Ns = 70; Ncirc = 20;
Ns = 25; Ncirc = 10;
%Ns = 300; Ncirc = 101;
%Ns = 230; Ncirc = 70;

dottype = 'antidot'; %dot or antidot (= '')
adivdL = 6; 
for adivd = adivdL
savepath = runLatticeHoneycomb('triangular','irrfbz',Ns,Ncirc,adivd,3,[],[],[],dottype);
continue
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

lnsty = {'.','-'};
lncol = {cols{14},cols{4}*.5+cols{8}*.5};

%plotting
try; close(2); end
set_figsize(2,3.5*numel(models)/2+3.5*.75,9*1.85*.75)

ylims = [0,.3536*models(1).ef_eV]*eV2THz;
yticks = 0.516958439800970*models(1).ef_eV*eV2THz*(0:4)/5;

subaxis(1,1,1,'SH',.01,'MR',.026,'ML',.165,'SV',0.014,'MT',.02,'MB',.08)
%plotting data
for nn = 2:numel(xticks)-1
    plot(xticks(nn)*[1,1],ylims,'-','color',cols{8},'linewidth',.4); hold on
end
for mm = 1 %2:-1:1
    %plot(kplot,real(eneeig_eV{mm}(:,shifti))*eV2THz,lnsty{mm},'color',lncol{mm},'linewidth',.8);
    plot(kplot,real(newene(:,shifti))*eV2THz,lnsty{mm},'color',lncol{mm},'linewidth',.8);
end
hold off

%annotations, clean up, etc
ylim(ylims)
xlim(minmax(kplot))
set(gca,'XTick',xticks,'XTickLabel',{'\bfM','\bf\Gamma','\bfK','\bfM'},'Fontsize',8,...
    'YTick',yticks,'LineWidth',.4,'TickLength',[.015,.015])

ylabel('Frequency (THz)','Fontsize',10)
xlabel('Wave vector  \bfk','Fontsize',10)

text(kplot(end)*.025,ylims(2)*.09,{['\itE\rm_F = ' num2str(models(mm).ef_eV) ' eV'],...
    ['\ita\rm = ' num2str(models(mm).L*1e9) ' nm']},...
    'Fontsize',8,'Color','k',...
    'VerticalAlignment','top','HorizontalAlignment','left','Margin',0.1)


%Inset of unit cell (Fine tuning by hand is necessary here)
handaxes2 = axes('Position', [0.55 -0 0.35 0.35]);

R1 = [1,0]; R2 = [cosd(60), sind(60)];
uc = [-R1-R2;  -R2+R1;  R2+R1; R2-R1]/2;
shift = [(1+cosd(60))/2,sind(60)/2]/(3);
theta = linspace(0,2*pi,51);

%Plot the unit cell
if strcmpi(dottype,'antidot') || strcmpi(dottype,'')
    bzcol = cols{4}*.15+cols{8}*.85;  dotcol = 'w';
elseif strcmpi(dottype(2:4),'dot')
    bzcol = 'w';                      dotcol = cols{4}*.15+cols{8}*.85;
    dottype = '_dot'; %prepare for saving
end
patch(uc([1:end,1],1),uc([1:end,1],2),bzcol,'linewidth',.75)
hold on
for pm = [-1,1]
    patch(cos(theta)/adivd/2+pm*shift(1),...
        sin(theta)/adivd/2+pm*shift(2),dotcol,'linewidth',.75)
end
set(handaxes2, 'box', 'off','color','none'); axis equal tight off
text(-0.65,-.6,['\itd\rm/\ita\rm = 1/' num2str((adivd),'%g')],'Fontsize',8);
drawnow


%save figure
%export_fig(['../figures/lattice/honeycomb' dottype '_adivd' strrep(num2str(adivd),'.','p')],'-pdf','-painters')
export_fig(['../figures/lattice/honeycomb' dottype '_adivd' strrep(num2str(adivd),'.','p')],'-pdf','-painters')
end