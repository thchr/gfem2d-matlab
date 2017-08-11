function DR = calcLatticeDifferentialAnisotropic(blochmesh,k)
%CALL:       DR = calcLatticeDifferentialAnisotropic(blochmesh,k)
%Calculates the anisotropic parts of the differential operator DR which 
%discretizes the FEM solution of the continuity equation on the mesh whose
%triangulation 'blochmesh' in a system of discrete translational symmetry 
%and with a k-vector ([2x1] column vector) in the FBZ.

%Unwrapping the necessary elements of the blochmesh structure
p = blochmesh.p;
t = blochmesh.t;
areas = blochmesh.area;
remesh = blochmesh.remesh;

if all(size(k) == [1,2]); %Ensure that k is a column vector
    k = k.';
end

%Check if k is a zero vector (zerok = 1 in that case, otherwise zerok = 0)
zerok = all(k == [0;0]);
R = [0,1;-1,0];

%Number of vertices and elements (triangles)
Nvert = size(p,1);
Ntri = size(t,1);

%Steps between each "progress" printout; we want 10 prints in total
printstep = floor(Ntri/10);


fprintf('Calculates the differential (anisotropic) matrix parts of DR')


%Construct right and left D matrices by looping
DR = sparse(Nvert,Nvert); 
timeD = tic;
for j = 1:Ntri
    if mod(j,printstep) == 0; fprintf('   j-loop:   %g/%g (%.1f min/%.1f min)\n',j,Ntri,toc(timeD)/60,toc(timeD)/60/(j/Ntri)); end
    
    drbc = remesh.p(remesh.t(j,2),:) - remesh.p(remesh.t(j,3),:); %Here we use the non-folded mesh, because it is the
    drca = remesh.p(remesh.t(j,3),:) - remesh.p(remesh.t(j,1),:); %direct distance (not the folded!) that enters in
    drab = remesh.p(remesh.t(j,1),:) - remesh.p(remesh.t(j,2),:); %the derivatives. Size is [1 x 2] row vectors
    
    dR = [ 0             , drbc*R*drca.' , drbc*R*drab.';
           drca*R*drbc.' , 0             , drca*R*drab.';
           drab*R*drbc.' , drab*R*drca.' , 0             ] / (4*areas(j));
        
        if ~zerok  %Add k-dependent parts if k is not zero
            dR2(1,1:3) = (drbc*k)*(1i/6);
            dR2(2,1:3) = (drca*k)*(1i/6);
            dR2(3,1:3) = (drab*k)*(1i/6);
                        
            dR = dR + dR2 + dR2.';
        end

    %Assigns right D components [in the blochmesh-basis!]
    DRj = sparse(Nvert,Nvert);   
    for n = 1:3
        for m = 1:3
            DRj(t(j,n),t(j,m)) = dR(n,m);
        end
    end
    
    %Sums contributions to total matrices
    DR = DR + DRj;
end