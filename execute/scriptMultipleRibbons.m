clc; clear all; close all;
addpath(genpath('..'))

%Some plotting options
cols=flatcolors;
plotopts = {':o','Color',cols{8},'MarkerFaceColor',cols{1},'MarkerEdgeColor',cols{8},'MarkerSize',5};

%Inclusion boundary
adivd = 2;
Ncirc = 30;
theta = linspace(0,2*pi,Ncirc+1); theta = theta(1:end-1);
inclusion = 1/(2*adivd)*[cos(theta);sin(theta)].';

%Calculation method for Coulomb matrix
CoulombMethod = '2Dsupercell'; %Can be '1Dseries' or '2Dsupercell'

%Specify setup parameters
setups{1} = struct('ncell',1,'Ns',70);
setups{2} = struct('ncell',2,'Ns',70);
setups{3} = struct('ncell',3,'Ns',100);
setups{4} = struct('ncell',4,'Ns',80);
setups{5} = struct('ncell',5,'Ns',70);
setups{6} = struct('ncell',6,'Ns',69);
setups{7} = struct('ncell',7,'Ns',61);
setups{8} = struct('ncell',8,'Ns',53);

for ss = 1:numel(setups)
    ncell = setups{ss}.ncell;
    Ns =    setups{ss}.Ns;
    
    %% MESHING
    blochmesh=geomRibbonInclusion('triangular', [1,ncell],inclusion,Ns,@(x,y) .015*250/Ns*ones(size(x)),'remesh');
    
    %Momentum list
    Nk = 25; %Number of k-points
    kvec = bsxfun(@times,linspace(0,1/2,Nk)',blochmesh.Grib); maxk = max(sqrt(kvec(:,1).^2 + kvec(:,2).^2));
    
    %% ACTUAL CALCULATION AND PLOTTING
    
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
        zeta(:,kk) = calcEigen(W,DL,DRi,numeigs,[],true);
        
        if plotfigs == 1
            %--- PLOTTING ON THE FLY ---
            figure(fighndl);
            plot(kvec(kk:Nk,1)/maxk,real(zeta(:,kk:Nk)),plotopts{:}) %Plot numerical computation
            xlim(minmax(kvec(:,1))/maxk); ylim([0,5/ncell])
            drawnow
        end
        
        %--- PRINT PROGRESS ---
        fprintf('DISPERSION LOOP %g/%g (%.2f/%.2f min)\n\n',Nk+1-kk,Nk,toc(tick)/60,toc(tick)/60*Nk/(Nk+1-kk));
    end
    %%
    if plotfigs == 1
        xlabel('\itk\rm/(\pi/\ita\rm)','Fontsize',11)
        ylabel('Dimensionless eigenvalue','Fontsize',11)
        set(gca,'Fontsize',8,'LineWidth',.1,'XTick',[0:.25:1])
        export_fig(['../figures/ribbons/ribbon_tri_ncell' num2str(ncell) '_' CoulombMethod '_Ns' num2str(Ns) '_Ncirc' num2str(Ncirc)],'-pdf')
    end
    
    clear W DRi DL
    save(['../output/ribbons/tri_ncell' num2str(ncell) '_' CoulombMethod '_Ns' num2str(Ns) '_Ncirc' num2str(Ncirc) '.mat'])
    %%
end