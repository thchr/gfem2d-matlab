function [DR,DL,D] = calcLatticeDifferential(blochmesh,k,f)
%CALL:          [DR,DL,D] = calcLatticeDifferential(blochmesh,k,f)
%Calculates the differential operators D, DL, and DR which discretize the
%FEM solution of the continuity equation on the mesh whose triangulation
%data are given by p (vertices) and t (triangle connections) in a system of
%discrete translational symmetry and with a k-vector ([2x1] column vector)
%in the first Brillouin zone.
%If f is supplied, the conductivity is assumed to be of an anisotropic form
%whose [xy] representation is
%       sigma = sigma_0*[ f(1,1), f(1,2); f(2,1), f(2,2) ],
%where sigma_0 is some chosen normalization such that f is dimensionless
%(with it being optimal to choose elements of f near unity).

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

%Number of vertices and elements (triangles)
Nvert = size(p,1);
Ntri = size(t,1);

%Steps between each "progress" printout; we want 10 prints in total
printstep = floor(Ntri/10);

%Isotropy/anisotropy setup
if ~exist('f','var') || (isnumeric(f) == 1 && f == 1) %Isotropic scenario
    anisotropy = 0;
    kflip = [-k(2);k(1)];
elseif all(size(f)==[2,2])               %Anisotropic scenario
    anisotropy = 1;
    frot1 = [ f(2,2) , -f(1,2) ; -f(2,1) , f(1,1) ]; %The "rotated" and transposed-like variants of the matrix f, which enter in the dR(a,c,d) elements
    frot3 = [ -f(2,1) , -f(2,2) ; f(1,1) , f(1,2) ]; frot3k = frot3*k;
    frot4 = [ -f(1,2) , -f(2,2) ; f(1,1) , f(2,1) ]; frot4k = frot4*k;
else                                     %Erroneous input
    error('Anisotropic calculations require f given as a 2x2 matrix')
end

fprintf('Calculates the differential matrices DR')
if nargout >= 2; fprintf(' and DL\n'); end; fprintf('\n')

%Construct right and left D matrices by looping
DR = sparse(Nvert,Nvert); DL = DR;
timeD = tic;
for j = 1:Ntri
    if mod(j,printstep) == 0; fprintf('   j-loop:   %g/%g (%.1f min/%.1f min)\n',j,Ntri,toc(timeD)/60,toc(timeD)/60/(j/Ntri)); end
    
    drbc = remesh.p(remesh.t(j,2),:) - remesh.p(remesh.t(j,3),:); %Here we use the non-folded mesh, because it is the
    drca = remesh.p(remesh.t(j,3),:) - remesh.p(remesh.t(j,1),:); %direct distance (not the folded!) that enters in
    drab = remesh.p(remesh.t(j,1),:) - remesh.p(remesh.t(j,2),:); %the derivatives. Size is [1 x 2] row vectors
    
    if anisotropy == 0; %Isotropic calculation
        dR = [ drbc*drbc.' , drbc*drca.' , drbc*drab.';
            drca*drbc.' , drca*drca.' , drca*drab.';
            drab*drbc.' , drab*drca.' , drab*drab.' ] / (4*areas(j));
        
        if ~zerok  %Add k-dependent parts if k is not zero
            dR2 = [2,1,1; 1,2,1; 1,1,2]*((k.'*k)*areas(j)/12);
            
            dR3(1:3,1) = (drbc*kflip)*(-1i/6);
            dR3(1:3,2) = (drca*kflip)*(-1i/6);
            dR3(1:3,3) = (drab*kflip)*(-1i/6);
            
            dR4(1,1:3) = (drbc*kflip)*(1i/6);
            dR4(2,1:3) = (drca*kflip)*(1i/6);
            dR4(3,1:3) = (drab*kflip)*(1i/6);
            
            dR = dR + dR2 + dR3 + dR4;
        end
    else %Anisotropic calculation
        dR = [ drbc*(frot1*drbc.') , drbc*(frot1*drca.') , drbc*(frot1*drab.');
            drca*(frot1*drbc.') , drca*(frot1*drca.') , drca*(frot1*drab.');
            drab*(frot1*drbc.') , drab*(frot1*drca.') , drab*(frot1*drab.') ] / (4*areas(j));
        
        if ~zerok  %Add k-dependent parts if k is not zero
            dR2 = [2,1,1; 1,2,1; 1,1,2]*((k.'*f*k)*areas(j)/12);
            
            dR3(1:3,1) = (drbc*frot3k)*(-1i/6);
            dR3(1:3,2) = (drca*frot3k)*(-1i/6);
            dR3(1:3,3) = (drab*frot3k)*(-1i/6);
            
            dR4(1,1:3) = (drbc*frot4k)*(1i/6);
            dR4(2,1:3) = (drca*frot4k)*(1i/6);
            dR4(3,1:3) = (drab*frot4k)*(1i/6);
            
            dR = dR + dR2 + dR3 + dR4;
        end
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
fprintf('\n')