function [zeta,eigV] = calcEigen(V,DL,DR,numeigs,method,keepzerosols)
%CALL   |   [zeta,eigV] = calcEigen(V,DL,DR,numeigs,method,keepzerosols)
%Calculate eigenvalues 'zeta' and vectors 'eigV' for the system described
%by V and DL and DR. Pass 'numeigs' as empty [] to retain ALL eigenvalues
%and vectors; otherwise enter a fixed number (output sorted from min-max)
%Possible 'method' choices are 'D' (fastest; default) and 'DLDRfull' (slow)

%Get version of Matlab (because eig seems to change behavior around v2013
vmatlab = version('-release'); vmatlab = vmatlab(1:4);

ticEigs = tic;
%Default values
if ~exist('numeigs','var')
    numeigs = 48;
end
if ~exist('method','var') || isempty(method)
    method = 'D';
end

%Solve eigensystem
switch method
    case 'D' %Compute D = DL\DR and solve a regular EVP rather than a generalized EVP (faster since DL and DR are sparse, so DL\DR is fast)
        if vmatlab > 2012
            [eigV,zeta]=eig(1/(2*pi)*(DL\DR)*V,'vector');  %Find ALL eigenvalues and vectors
        else
            [eigV,zeta]=eig(1/(2*pi)*(DL\DR)*V); zeta = diag(zeta);
        end
    case 'DLDRfull' %Solve as generalized EVP (full, unfortunately)
        [eigV,zeta]=eig(1/(2*pi)*DR*V,full(DL),'vector');  %Find ALL eigenvalues and vectors
end
% [eigV,ZETA]=eigs(1/(2*pi)*DR*V,DL,numeigs,'sr');  %Find numeigs with smallest possible real part

%Sort according (small to large)
[zeta,I]=sort(zeta); eigV = eigV(:,I);

%Remove any spurious near-zero eigenvalue solutions (arising due to discretization)
if ~exist('keepzerosols','var') || isempty(keepzerosols) || keepzerosols == 0 %Default choice is to remove; if keepzerosols = 1 we do NOT
    eigV(:,real(zeta)<1e-7)  = []; zeta(real(zeta)<1e-7)  = [];               %remove anything (e.g. meaningful at Gamma point in lattices)
end

%Retain only the numeig lowest eigenstates
if ~isempty(numeigs)
    zeta = zeta(1:numeigs); eigV = eigV(:,1:numeigs);
end

%Print how long it took
fprintf('   Eigenvalues/vectors computed in %.1f min (%g lowest eigenpartners output)\n',toc(ticEigs)/60, numel(zeta))

