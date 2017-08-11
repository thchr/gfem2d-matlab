clear all; close all; clc;
ConstantsUnits0;
ef_eV = .2; L = 100e-9;

gam_eV = 6e-3;
BL = [0,5,10];
omegac_eV = e*vf^2*BL/(ef_eV*ev2jo)*hbar_eV;

%models(1) = struct('type','isotropic','ef_eV',ef_eV,'L',L,'B',BL(1),'system',[]        ,'approx','intra','keep_eV',[omegac_eV(1),ef_eV],'gam_eV',gam_eV,'gamc',[]);
%models(2) = struct('type','magneto'  ,'ef_eV',ef_eV,'L',L,'B',BL(2),'system','graphene','approx',[]     ,'keep_eV',[omegac_eV(2),ef_eV],'gam_eV',[]    ,'gamc',gam_eV/omegac_eV(2)); %The lower cutoff is the cyclotron frequency
%models(3) = struct('type','magneto'  ,'ef_eV',ef_eV,'L',L,'B',BL(3),'system','graphene','approx',[]     ,'keep_eV',[omegac_eV(3),ef_eV],'gam_eV',[]    ,'gamc',gam_eV/omegac_eV(3));
%savepath = runRibbonLatticeMagneto([],10,41,30,2,49,0,models);

savepath = '../output/ribbons/magnet_triangular_numcell10_Ns41_Ncirc30_savephi0.mat';
%%
load(savepath)
eneeig_eV{1}(imag(eneeig_eV{1})>0) = NaN;
eneeig_eV{2}(real(eneeig_eV{2})<omegac_eV(2)*1.01) = NaN;
eneeig_eV{3}(real(eneeig_eV{3})<omegac_eV(3)*1.01) = NaN;
%%
ene_eV = linspace(0,.965*ef_eV,5000);
for mm = 1:3; dos{mm} = calcDOS(eneeig_eV{mm},ene_eV,L); end
maxdos = max( [ max(dos{1}(:)), max(dos{2}(:)), max(dos{3}(:)) ]);

%%

cols=flatcolors;

set_figsize(1,18,13)
cmap = bsxfun(@minus,[1,1,1],flipud(colormap(gray))); %Setting up a white-to-cols{*} colormap
cmap = bsxfun(@minus,[1,1,1],bsxfun(@times,cmap,[1,1,1]-cols{14}));
for mm = 1:3
    
    
    subaxis(1,3,mm,'SH',.015,'MR',.01,'ML',.075,'MT',.05,'MB',.08)
    imagesc(kvec(:,1)/maxk,ene_eV,dos{mm}/1e9);
    %plot(kvec(:,1)/maxk,real(eneeig_eV{mm}),'.','color',cols{8},'MarkerSize',2)
    %contourf(kvec(:,1)/maxk,ene_eV,dos{mm}/1e9,linspace(0,maxdos/1e9,50),'LineColor','none');
    %shading interp
    caxis([0,maxdos/1e9*.70])
    
    colormap(cmap)
    %colormap(morgenstemning)
    if mm ~= 1
        patch([-1,1,1,-1],[omegac_eV(mm),omegac_eV(mm),0,0],cols{4}*.15+[1,1,1]*.85,'LineStyle','none')
    else
        cb=colorbar([0.3150396475770924 0.1 0.0102790014684288 0.124489795918367],'FontSize',6,'YTick',[0,5,10,15]);
        text(.525,.019,'DOS [eV^{-1}nm^{-1}]','Rotation',90,'HorizontalAlignment','center','Fontsize',7,'LineWidth',.1)
    end
    set(gca,'Layer','Top','Fontsize',8,'LineWidth',.1,'YDir','normal');
    if mm ~= 1; set(gca,'YTickLabel',{''}); else ylabel('Energy [eV]','Fontsize',10);  end
    if mm ==2; xlabel('Momentum \itk_x\rm/(\pi/\ita\rm)','Fontsize',10); end
    xlim([-1,1])
    ylim([minmax(ene_eV)])
    title(['\itB\rm = ' num2str(BL(mm)) ' T'],'Fontsize',10)
    drawnow; pause(0.01);
end

%export_fig('DOS_ncell10_ef200meV_gam6meV_a100nm_adivd2_pix','-pdf')