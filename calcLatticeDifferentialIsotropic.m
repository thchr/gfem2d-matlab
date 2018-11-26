function [DR,DL,D] = calcLatticeDifferentialIsotropic(blochmesh,k)
%CALL:       [DRani] = calcLatticeDifferentialIsotropic(blochmesh,k)
%Calculates the isotropic parts of the differential operator DR, DL, and D  
%which discretizes the FEM solution of the continuity equation on the mesh 
%whosetriangulation 'blochmesh' in a system of discrete translational  
%symmetry and with a k-vector ([2x1] column vector) in the FBZ.

%Unwrapping the necessary elements of the blochmesh structure
p = blochmesh.p;
t = blochmesh.t;
areas = blochmesh.area;
remesh = blochmesh.remesh;

if all(size(k) == [1,2]) %Ensure that k is a column vector
    k = k.';
end

%Check if k is a zero vector (zerok = 1 in that case, otherwise zerok = 0)
zerok = all(k == [0;0]);
kflip = [-k(2);k(1)];

%Number of vertices and elements (triangles)
Nvert = size(p,1);
Ntri = size(t,1);

%Steps between each "progress" printout; we want 10 prints in total
printstep = floor(Ntri/10);

fprintf('Calculates the differential (isotropic) matrix parts of DR')
if nargout >= 2; fprintf(' and DL'); end; fprintf('\n')

%Construct right and left D matrices by looping
DR = sparse(Nvert,Nvert); DL = DR;
timeD = tic;
for j = 1:Ntri
    if mod(j,printstep) == 0; fprintf('   j-loop:   %g/%g (%.1f min/%.1f min)\n',j,Ntri,toc(timeD)/60,toc(timeD)/60/(j/Ntri)); end
    
    drbc = remesh.p(remesh.t(j,2),:) - remesh.p(remesh.t(j,3),:); %Here we use the non-folded mesh, because it is the
    drca = remesh.p(remesh.t(j,3),:) - remesh.p(remesh.t(j,1),:); %direct distance (not the folded!) that enters in
    drab = remesh.p(remesh.t(j,1),:) - remesh.p(remesh.t(j,2),:); %the derivatives. Size is [1 x 2] row vectors
    
    dR = [ drbc*drbc.' , drbc*drca.' , drbc*drab.';
           drca*drbc.' , drca*drca.' , drca*drab.';
           drab*drbc.' , drab*drca.' , drab*drab.' ] / (4*areas(j));
        
        if ~zerok  %Add k-dependent parts if k is not zero
            dR2(1,1:3) = (drbc*kflip)*(1i/6);
            dR2(2,1:3) = (drca*kflip)*(1i/6);
            dR2(3,1:3) = (drab*kflip)*(1i/6);
            
            dR3 = [2,1,1; 1,2,1; 1,1,2]*((k.'*k)*areas(j)/12);
            
            dR = dR + dR2 - dR2.' + dR3;
        end

    
    
    %Calculate left and right D 3x3 blocks
    if nargout >= 2;  dL = [2,1,1; 1,2,1; 1,1,2]*areas(j)/12; end
    
    %Assigns left and right D components [in the blochmesh-basis!]
    DRj = sparse(Nvert,Nvert);   if nargout >= 2; DLj = sparse(Nvert,Nvert); end
    for n = 1:3
        for m = 1:3
            DRj(t(j,n),t(j,m)) = dR(n,m);
            if nargout >= 2; DLj(t(j,n),t(j,m)) = dL(n,m); end
        end
    end
    
    %Sums contributions to total matrices
    DR = DR + DRj;
    if nargout >= 2; DL = DL + DLj; end
end

%Calculate total D matrix by inversion
if nargout >= 3;
    fprintf('Calculating the full D-matrix by solving -D_L\\D_R\n')
    timeInv=tic;
    D = -full(DL\DR);
    fprintf('   Inversion completed in %.1f min\n',toc(timeInv)/60)
end