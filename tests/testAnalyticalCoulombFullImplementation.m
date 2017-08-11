clc; clear all; close all;

addpath(genpath('..'))
%Mesh settings

zetaexact = [1.0977 4.9140 8.1337 11.3079;
    1.9942 6.2455 9.5455 12.7592;
    2.8556 7.5124 10.8989 14.1596;
    3.7032 8.7395 12.2117 15.5221];
zetaexact = reshape(zetaexact,16,1);



edgedens = 50:50:500;
for mm=1:numel(edgedens)
    ticmesh=tic;
    [~,p,t] = evalc('geomDisk(edgedens(mm),15/edgedens(mm),1)');
    fprintf('Mesh with %g nodes and %g elements (created in %.1f sec)\n',size(p,1),size(t,1),toc(ticmesh))
    meshval(mm) = size(p,1);
    %%
    
    t0 = tic; 
    [~,Vn] = evalc('calcCoulomb(p,t,''numquad'')');    
    [~,Va] = evalc('calcCoulomb(p,t,''analytical'')');
    [~,DR,DL] = evalc('calcDifferential(p,t)');
    fprintf('   V, DL, and DR calculated in %.1f sec\n',toc(t0))
    
    fprintf('   numquad    |')
    zetan{mm} = calcEigen(Vn,DL,DR);
    
    fprintf('   analytical |')
    zetaa{mm} = calcEigen(Va,DL,DR);
    
    fprintf('   Mesh %g/%g CALCULATED (total time: %.1f sec)\n\n\n',mm,numel(edgedens),toc(ticmesh))
end
%%

N = 35;
[cols,nams] = flatcolors;
set_figsize([],17,17)
hold on
for zz = 1:numel(zetaexact)
    pE=plot([0,6000],zetaexact(zz)*[1,1],'-','color',cols{4});
end
for mm = 1:numel(meshval)
    pQ=plot(meshval(mm)*ones(1,N),zetan{mm}(1:N),'.','Color',cols{14}*.6+cols{24}*.4,'Markersize',8.5);
    pA=plot(meshval(mm)*ones(1,N),zetaa{mm}(1:N),'o','Color',cols{1},'MarkerSize',3,'LineWidth',.95);
end
hold off
ylim([0,7.75])
legend([pE,pQ,pA],{'Exact \itl\rm \neq 0 eigenvalues','''numquad'' Coulomb','''analytical'' Coulomb'},'Fontsize',9,'Location','SouthEast'); legend boxoff
xlabel('Number of nodes','Fontsize',10)
ylabel('Eigenvalues \zeta','Fontsize',10)
box on
set(gca,'Fontsize',7,'LineWidth',.125)

export_fig('ConvergenceDisk_NodesVsZeta','-pdf')