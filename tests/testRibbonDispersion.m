clc; clear all; close all;
cols=flatcolors;
%addpath(genpath('..'))

%Inclusion boundary
adivd = 2;
Ncirc = 30;
theta = linspace(0,2*pi,Ncirc+1); theta = theta(1:end-1);
inclusion = 1/(2*adivd)*[cos(theta);sin(theta)].';


ncell = 8;
Ns = 53;
%blochmesh=geomRibbonInclusion({[1,0],[0,1]},[1,1],[],Ns,@(x,y) .015*250/Ns*ones(size(x)),'remesh');
blochmesh=geomRibbonInclusion('triangular', [1,ncell],inclusion,Ns,@(x,y) .015*250/Ns*ones(size(x)),'remesh');

%Momentum list
Nk = 25; %Number of k-points
kvec = bsxfun(@times,linspace(0,1/2,Nk)',blochmesh.Grib); maxk = max(sqrt(kvec(:,1).^2 + kvec(:,2).^2));
%kvec = blochmesh.Grib/2;
CoulombMethod = '1Dseries'; %Can be '1Dseries' or '2Dsupercell'

plotopts = {':o','Color',cols{8},'MarkerFaceColor',cols{1},'MarkerEdgeColor',cols{8},'MarkerSize',5};
%%
c(1,:) = [.8216, 0     , -.6603, .4648, .5982];
c(2,:) = [.8216, 2.3159, -.1787, .2761, .1334];
c(3,:) = [1    , 5.5103,  .0378, .1723, .0376];
c(4,:) = [1    , 8.6348,  .1787, .1305, .0144];
zetaa = @(mm,x) c(mm,1)*x + (c(mm,2) + c(mm,3)*x )./(1 + c(mm,4)*x + c(mm,5)*x.^2);

try close(2); end
fighndl= set_figsize(2,18,18);

numeigs = 48;
tick = tic;
for kk = Nk:-1:1
    %--- ACTUAL CALCULATION ---
    switch CoulombMethod
        case '1Dseries'
            W = calcRibbonCoulombViaSeries(blochmesh,kvec(kk,:),10);
        case '2Dsupercell'
            if kk == Nk; blochmesh.R{2} = 3*ncell*blochmesh.R{2}; end
            W = calcLatticeCoulombOptimized(blochmesh,kvec(kk,:));
    end
    
    [DRi,DL] = calcLatticeDifferentialIsotropic(blochmesh,kvec(kk,:));
    [zeta(:,kk),eigV(:,:,kk)] = calcEigen(W,DL,DRi,numeigs,[],true);
    
    %--- PLOTTING ON THE FLY ---
    figure(fighndl);
  %  plot(3/4*[1,1],[0,1000],':','Color',cols{4}); hold on;
    plot(kvec(kk:Nk,1)/maxk,real(zeta(:,kk:Nk)),plotopts{:}) %Numerical computation
    
    hold on
    %Semi-analytical from fit
    for mm = 1:size(c,1) %Different modes in radial/width sense
        for ff = 0:4     %Zone folding
            if mod(ff,2) == 0; x = kvec(:,1)/maxk; else x = flipud(kvec(:,1))/maxk; end
            plot(x,zetaa(mm,kvec(:,1)+ff*norm(blochmesh.Grib,2)/2),'-','color',cols{mm+1})
        end
    end; 
    
    hold off
    xlim(minmax(kvec(:,1))/maxk); ylim([0,5/ncell])
    drawnow
    
    
    %--- PRINT PROGRESS ---
    fprintf('DISPERSION LOOP %g/%g (%.2f/%.2f min)\n\n',Nk+1-kk,Nk,toc(tick)/60,toc(tick)/60*Nk/(Nk+1-kk));
end
%%
%export_fig('../figures/ribbon_1Dvia2Dsuperlattice_check','-pdf')
%export_fig('../figures/ribbon_1D_check','-pdf')
xlabel('\itk\rm/(2\pi/\ita\rm)','Fontsize',11)
ylabel('Dimensionless eigenvalue','Fontsize',11)
set(gca,'Fontsize',8,'LineWidth',.1,'XTick',[0:.25:1])
export_fig(['../figures/ribbon_tri_ncell' num2str(ncell) '_' CoulombMethod '_Ns' num2str(Ns) '_Ncirc' num2str(Ncirc)],'-pdf')
save(['../output/ribbons/tri_ncell' num2str(ncell) '_' CoulombMethod '_Ns' num2str(Ns) '_Ncirc' num2str(Ncirc) '.mat'])
%%
break

%%
p = blochmesh.remesh.p;
t = blochmesh.remesh.t;
try; close(3); end
set_figsize(3,40,20);
for nn = 1:12
    subaxis(4,3,nn,'M',.01,'MT',.05,'SV',0.04,'SH',0.02)
    
    v = blochmesh.remesh.values(W*eigV(:,nn,kk));
    for rr=-1:1;
        %plotVertexData(p,t,imag(v));
        pR = bsxfun(@plus,p,rr*blochmesh.Rrib);
        vtot = v.*exp(1i*kvec(kk,:)*pR.').';
        trisurf(t,pR(:,1),pR(:,2),real(vtot),'EdgeColor','none'); hold on
    end
    hold off
    shading interp;
    colormap('bluewhitered_mod'); freezeColors;
    view(2);
    axis equal off;
    drawnow;
    title(zeta(nn,kk))
end