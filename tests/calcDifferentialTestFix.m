function [DR,DL,D] = calcDifferentialTestFix(p,t)
%CALL:          [DR,DL,D] = calcDifferentialTestFix(p,t)
%DESCRIPTION: Calculates the differential operators D, DL, and DR which
%discretize the FEM solution of the continuity equation on the mesh whose
%triangulation data are given by p (vertices) and t (triangle connections).
%If f is supplied, the conductivity is assumed to be of an anisotropic form
%whose [xy] representation is
%       sigma = sigma_0*[ f(1,1), f(1,2); f(2,1), f(2,2) ],
%where sigma_0 is some chosen normalization such that f is dimensionless
%(with it being optimal to choose elements of f near unity).


%Mesh info
Nvert = size(p,1);
Ntri = size(t,1);
area = meshArea(p,t);

fprintf('Calculates the differential matrices DR and DL (if req''ed)\n')


%Construct right and left D matrices by looping
timeD = tic;
D = sparse(Nvert,Nvert); DL = D; DR = D; 
for j = 1:Ntri
    if mod(j,500) == 0; fprintf('   j-loop:   %g/%g (%.1f min/%.1f min)\n',j,Ntri,toc(timeD)/60,toc(timeD)/60/(j/Ntri)); end
    
    
    drbc = p(t(j,2),:) - p(t(j,3),:); drca = p(t(j,3),:) - p(t(j,1),:); drab = p(t(j,1),:) - p(t(j,2),:);
    %Isotropic calculation
    tn = t; tn(j,:) = NaN; %([1:tt-1,tt+1:end],:);
    
    dR = zeros(3,3); ec = 1; 
    if any(all([any(tn == t(j,2),2), any(tn == t(j,3),2)],2)) %side 1
        dR = dR + [ 0 , drbc*drbc.' , drbc*drbc.' ; ...
                    0 , drca*drbc.' , drca*drbc.' ; ...
                    0 , drab*drbc.' , drab*drbc.'  ].';
    else
        ec = 0; 
    end
    if any(all([any(tn == t(j,3),2), any(tn == t(j,1),2)],2)) %side 2
        dR = dR + [ drbc*drca.' , 0 , drbc*drca.' ; ...
                    drca*drca.' , 0 , drca*drca.' ; ...
                    drab*drca.' , 0 , drab*drca.'  ].';
    else
        ec = 0; 
    end
    if any(all([any(tn == t(j,1),2), any(tn == t(j,2),2)],2)) %side 3
        dR = dR + [ drbc*drab.' , drbc*drab.' , 0 ; ...
                    drca*drab.' , drca*drab.' , 0 ; ...
                    drab*drab.' , drab*drab.' , 0  ].';
    else
        ec = 0; 
    end
    if ec == 0;     
        dR
    end
    dR = -dR/(4*area(j));
    
    %         dR = [ drbc*drbc.' , drbc*drca.' , drbc*drab.';
    %                drca*drbc.' , drca*drca.' , drca*drab.';
    %                drab*drbc.' , drab*drca.' , drab*drab.' ] /(4*area(j));
    
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