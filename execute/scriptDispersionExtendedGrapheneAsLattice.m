clear all; clc; close all;

addpath(genpath('..'));
ConstantsUnits0; [cols,nams]=flatcolors;

%[p,t] = geomDisk();
%[p,t] = geomAnnulus();drawnow;
%[p,t] = geomEquiTriangle(1,.04);drawnow;
%[p,t] = geomHoleSquareArray(4,1,[],[],struct('fun',@(x,y) .0125.*(sqrt(x.^2 + y.^2).^1.5+2)));
%[p,t] = geomHoleSquareArray(3,1);
%[p,t] = geomSquare(3,3,.04);
%blochmesh = geomPeriodicInclusion('square','disk',200,'antidot');

Ne = 100;
blochmesh = geomPeriodicInclusion('square',[],Ne,'const');
R = blochmesh.R;
G = calcReciprocal(R);

ef_eV = 0.2;
gam_eV = (1/6.37e-12)*hbar_eV;
T = 0;


%Momentum list
[kvec,kplot,kmark] = irreducibleFBZ('square',30,'plotfbz');

ticEL=tic;
for kk = 1:size(kvec,1);
    [zeta(:,kk)] = calcLatticeEigen(blochmesh,kvec(kk,:));
    fprintf('\n\n     DISPERSION LOOP: %g/%g (%.1f min/%.1f min)\n\n\n',kk,size(kvec,1),toc(ticEL)/60,size(kvec,1)/kk*toc(ticEL)/60);
end


%%
set_figsize([],14,14); 

plot(kplot, sqrt(zeta),'-','color',cols{2},'LineWidth',1.5)
hold on
for n = -5:5
    for m = -5:5
        Gnm = G{1}*n + G{2}*m;       
        val = sqrt( (kvec(:,1) + Gnm(1)).^2 + (kvec(:,2) + Gnm(2)).^2);
        plot(kplot,sqrt(val),':','color',cols{4});
    end
end
for mm = 2:size(kmark.n,1)-1
    plot(kplot(kmark.n(mm,1))*[1,1],[0,1]*7,':k')
end
hold off

box on
set(gca,'Layer','Top','Fontsize',8)
set(gca,'XTick',[kplot(kmark.n(:,1))],'XTickLabel',{kmark.symbol{kmark.n(:,2)}})
axis tight
ylim([0,4])
xlim([0,max(kplot)])
xlabel('Momentum','Fontsize',10)
ylabel('Dimensionless resonance value','Fontsize',10)
export_fig(['../figures/DispersionExtendedGraphene_Ne' num2str(Ne)],'-pdf')