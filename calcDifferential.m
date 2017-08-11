function [DR,DL,D] = calcDifferential(p,t,f)
%CALL:          [DR,DL,D] = calcDifferential(p,t,f)
%DESCRIPTION: Calculates the differential operators D, DL, and DR which 
%discretize the FEM solution of the continuity equation on the mesh whose
%triangulation data are given by p (vertices) and t (triangle connections). 
%If f is supplied, the conductivity is assumed to be of an anisotropic form
%whose [xy] representation is 
%       sigma = sigma_0*[ f(1,1), f(1,2); f(2,1), f(2,2) ],
%where sigma_0 is some chosen normalization such that f is dimensionless
%(with it being optimal to choose elements of f near unity).

fprintf('Calculates the differential matrices DR and DL (if req''ed)\n')
%Mesh info
Nvert = size(p,1);
Ntri = size(t,1);
area = meshArea(p,t);

%Isotropy/anisotropy setup
if ~exist('f','var') || (isscalar(f) && f==1)     %Isotropic scenario
    anisotropy = 0; 
elseif ismatrix(f) && all(size(f)==[2,2])         %Anisotropic scenario
    anisotropy = 1;
    frot = [ f(2,2) , -f(2,1) ; -f(1,2) , f(1,1) ]; %The "rotated" and transposed-like variant of the matrix f, which enter in the dR elements
    fprintf('   (input ''f'' is a 2x2 matrix, indicating anisotropy: f = [%g,%g;%g,%g])\n',f(1,1),f(1,2),f(2,1),f(2,2))
elseif ndims(f) == 3 
    fprintf('   |(input ''f'' is a 3-dimensional array: treated as spatial dependence,\n   |i.e. inhomogeneous material properties)\n')
    anisotropy = 2;                               %Anisotropic, inhomogenous case
    frott(1,1,:) = sum(f(2,2,t),2);   frott(1,2,:) = -sum(f(2,1,t),2); 
    frott(2,1,:) = -sum(f(1,2,t),2);  frott(2,2,:) = sum(f(1,1,t),2); 
    frott = frott/3;
else                                              %Erroneous input
    error('Anisotropic calculations require f given as a 2x2 matrix')
end

%Steps between each "progress" printout; we want 10 prints in total
printstep = floor(Ntri/10);

%Construct right and left D matrices by looping
timeD = tic; 
DR = sparse(Nvert,Nvert); if nargout>=2; DL = DR; end %Preallocate
for j = 1:Ntri
    if mod(j,printstep) == 0; fprintf('   j-loop:   %g/%g (%.1f min/%.1f min)\n',j,Ntri,toc(timeD)/60,toc(timeD)/60/(j/Ntri)); end

    
    drbc = p(t(j,2),:) - p(t(j,3),:); drca = p(t(j,3),:) - p(t(j,1),:); drab = p(t(j,1),:) - p(t(j,2),:);
    if anisotropy == 0; %Isotropic calculation
        dR = [ drbc*drbc.' , drbc*drca.' , drbc*drab.'; 
               drca*drbc.' , drca*drca.' , drca*drab.'; 
               drab*drbc.' , drab*drca.' , drab*drab.' ] /(4*area(j));
    else                %Anisotropic calculation
        if anisotropy == 2; 
            frot = frott(:,:,j); %Anisotropic, inhomogeneous calculation
        end
        dR = [ drbc*(frot*drbc.') , drbc*(frot*drca.') , drbc*(frot*drab.'); 
               drca*(frot*drbc.') , drca*(frot*drca.') , drca*(frot*drab.'); 
               drab*(frot*drbc.') , drab*(frot*drca.') , drab*(frot*drab.') ] /(4*area(j));
    end

    
    %Calculate left and right D 3x3 blocks
    if nargout >= 2;  dL = [2,1,1; 1,2,1; 1,1,2]*area(j)/12; end
    
    %Assigns left and right D components
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