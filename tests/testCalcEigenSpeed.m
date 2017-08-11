clc; clear all; close all;

addpath(genpath('..'))
%Mesh settings

edgedens = 40:20:300;
for mm=1:numel(edgedens)
    [~,p,t] = evalc('geomDisk(edgedens(mm),15/edgedens(mm))');
    meshval(mm) = size(p,1);
    fprintf('Mesh with %g nodes and %g elements\n',size(p,1),size(t,1))
    %%
    
    [~,V] = evalc('calcCoulomb(p,t,''analytical'')');
    [~,DR,DL] = evalc('calcDifferential(p,t)');
    
    t1=tic;
    zetaa{mm} = calcEigen(V,DL,DR,48,'DLDRfull');
    fprintf('     analytical Coulomb eigenvalues computed in %.1f sec\n',toc(t1))
    t2=tic;
    zetaa{mm} = calcEigen(V,DL,DR,48,'Dfull');
    fprintf('     analytical Coulomb eigenvalues computed in %.1f sec\n',toc(t2))
    
    fprintf('     Mesh %g/%g CALCULATED\n\n\n',mm,numel(edgedens))
end