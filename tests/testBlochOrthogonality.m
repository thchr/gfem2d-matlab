clear; close all; clc;


Ncirc = 36; Ns = 70; adivd=1.9; latticetype='triangular';
theta = linspace(0,2*pi,Ncirc+1); theta = theta(1:end-1); %Inclusion boundary
inclusion = 1/(2*adivd)*[cos(theta);sin(theta)].';        %----||----
mesh = geomPeriodicInclusion(latticetype,inclusion,Ns,@(x,y) .015*250/Ns*ones(size(x)),1);

bloch.type = '2D'; bloch.k = [0,.2];

%

%%
[~,mesh.remesh.areavec] = meshArea(mesh.remesh.p,mesh.remesh.t);
[~,V,DL,DRi] = evalc('calcSystemMatrices(mesh,bloch)');

%%
%clc
% solve eigenvalue problem
[~,ene,rho] = evalc('calcEigen(V,DL,DRi,5,[],all(bloch.k==0))'); % retain zero-solution if k = [0,0]

% get dual eigenvectors (potentials)
phi = V*rho;


% normalize the direct and dual eigenvectors st. <phi_n|rho_m> = delta_nm
for mm = 1:size(rho,2)
    overlap(mm,1) = integrateMeshFunction(mesh.remesh.p,mesh.remesh.t,...
        conj(mesh.remesh.values(rho(:,mm))).*mesh.remesh.values(phi(:,mm)),...
        mesh.remesh.areavec);
end
rho = bsxfun(@times,rho,1./sqrt(overlap).');
phi = bsxfun(@times,phi,1./sqrt(overlap).');

%Check that the states are orthonormal now
for mm = 1:size(rho,2)
    for nn = 1:size(rho,2)
        T(mm,nn) = integrateMeshFunction(mesh.remesh.p,mesh.remesh.t,...
            conj(mesh.remesh.values(rho(:,mm))).*mesh.remesh.values(phi(:,nn)),...
            mesh.remesh.areavec);
    end
end

disp(abs(T))

