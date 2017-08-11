function gradV = calcGradCoulombGeneralPoints(p,t,r)
%CALL:        gradV = calcGradCoulombGeneralPoints(p,t,r)
%DESCRIPTION: Calculates the gradCoulomb matrix for interaction between a 
%triangular mesh (p,t) and arbitrary (in-plane) points r (which is an Nx2
%array for xy-coords or an Nx3 array for xyz coords). Multiplication of
%gradV onto rho gives the electric field at point(s) r.
%IMPLEMENTATION: Uses the Wilton scheme for maximal accuracy. See:
%   'WORK NOTE: Coulomb integrals over triangular 2D mesh elements:
%   applications to FEM-BEM method for graphene plasmonics' 
%   Wilton, Rao, Glisson, Schaubert, Al-Bundak, and Butler, 
%   'Potential integrals for uniform and linear source distributions on 
%   polygonal and polyhedral domains', Section VI.
%Created Dec. 21, 2016.

%Number of vertices and elements (triangles)
Nvert = size(p,1);
Ntri = size(t,1);

fprintf('Calculates the Coulomb matrix by analytical evaluation of integral terms (Wilton-scheme) V\n')

%Preallocation
gradV = zeros(size(r,1),Nvert);

timeV = tic;
for j = 1:Ntri
    if mod(j,500) == 0; fprintf('   j-loop:   %g/%g (%.1f min/%.1f min)\n',j,Ntri,toc(timeV)/60,toc(timeV)/60/(j/Ntri)); end
    l = t(j,1); m = t(j,2); n = t(j,3);
    rl = p(l,:); rm = p(m,:); rn = p(n,:);
    rtri = [rl;rm;rn]; %Vertices of j'th triangle (coordinates as rows [x y; x y; xy])
    
    %Each element is calculated by numerical evaluation of the analytical
    %expressions obtained by Wilton et al., and then added into the total
    %Coulomb matrix at the vertex sites
    gradV(:,[l m n]) = gradV(:,[l m n]) + real(calcGradWiltonContribution(rtri,r(:,1:2),r(:,3))); 
    %Small rounding errors in the logarithms (I suspect) may give rise to small
    %imaginary parts (near the machine error): we purge them immediately.
    %(in addition, this keeps V as a real double, which is advantageous for
    %addition - otherwise, the procedure may be significantly slowed)
end


fprintf('\n')

function gradVtri = calcGradWiltonContribution(rtri,r,z)
%The scheme employed here is explicated in the work note 'Coulomb integrals
%over triangular 2D mesh elements: application to FEM-BEM method for
%graphene plasmonics'. rtri is expected in the form [rl;rm;rn] = [r1;r2;r3]

%Construct the triangle-specific quantities peta & pxi
Rt = [0,-1;1,0];
denom = ((rtri(3,:) - rtri(1,:))*Rt)*(rtri(1,:)-rtri(2,:)).';
peta = -((rtri(3,:)-rtri(1,:))*Rt)/denom;
pxi = -((rtri(1,:)-rtri(2,:))*Rt)/denom;

%Construct the (transposed) a-vector and b-scalar (the latter, being r-dependent, is evaluated for at r-point, and the result stored in distinct columns)
r1mr = bsxfun(@minus,rtri(1,:),r).';
a(1,:) = -(peta+pxi);   b(1,:) = 1 - a(1,:)*r1mr;    
a(2,:) = peta;          b(2,:) = -a(2,:)*r1mr;
a(3,:) = pxi;           b(3,:) = -a(3,:)*r1mr;

%Construct the auxiliary parameters for evaluation of I0 and I1
rpm = [circshift(rtri,-1).', rtri.']; %rpm = [rminus;rplus], i.e. rpm(:,1:3) = r^+, rpm(:,4:6) = r^-. This ordering follows for most 'pm' postscripted quantities

l = bsxfun(@times,(rpm(:,1:3) - rpm(:,4:6)),1./sqrt(sum((rpm(:,1:3)-rpm(:,4:6)).^2,1))); l = [l,l];%Edge tangential unit vector
u = [l(2,:);-l(1,:)]; %Edge normal unit vector
alpha = [1,1,1,-1,-1,-1]; %Plus/minus prefactors 

I0 = zeros(size(r,1),1); gradI1 = zeros(size(r,1),3,2); gradI0 = zeros(size(r,1),3);
for ee = 1:6 %This part of the implementation could probably be vectorized and hence improved; too much effort at the moment
    drrpm = [bsxfun(@minus, r.', rpm(:,ee))]; %r - rpm
    lpm = -drrpm.'*l(:,ee); 
        
    P0u = drrpm.'*u(:,ee); %I use the oppposite sign convention for P0u compared to the notes here (and for the items it acts on) 
    R02 = P0u.^2 + z.^2; R0 = sqrt(R02);  
    Rpm = sqrt( sum( rmrplus.^2, 1 )  + z.^2.' ).';

    %Direct terms
    atanterm = atan(abs(P0u).*lpm ./(R02 + abs(z).*Rpm ));
    zterm = abs(z)./abs(P0u) .* atanterm;
    logterm = log( Rpm + lpm );
    logterm( isnan(logterm) | isinf(logterm) ) = 0; %Remove singular terms which are later multiplied by zero
    
    %Adding into the integral term
    I0 = I0 - alpha(ee) .* P0u .* ( logterm - zterm ); 
    
    
    sgnP0u = P0u./abs(P0u);
    for rr = 1:size(r,1); %For now, we loop over the input coordinates for simplicity
    %Gradient terms 
    gradR0 = ([P0u(rr)*u(:,ee);z(rr)])./R0(rr); 
    gradP0 = sgnP0u(rr)*[u(:,ee);0];
    gradRpm = ([drrpm(:,rr); z(rr)])./Rpm(rr);

    %Fraction term
    F = ( -P0.*lpm .* ( 2*R0(rr)*gradR0 + Rpm*sgn(z(rr))*[0;0;1] + abs(z(ee))*gradRpm ) + ...
          ( R02(rr) + abs(z)*Rpm(rr) ) * ( lpm(rr)*gradP0 - abs(P0u(rr))*l(:,ee) ) ) / ...
        ( (P0u(rr)*lpm(rr))^2 + ( R02(rr) + abs(z(rr))*Rpm(rr) )^2 ); 
    
    %Adding into the gradient integral terms
    gradI1 = gradI0 + alpha(ee) .* ( NaN );  %Not yet done !!!!!!!!!!!!!!!!!!!!!
    gradI0 = gradI0 + alpha(ee) .* ( ...
                - logterm(rr)*u(:,ee) ...
                + P0u(rr)/(Rpm(rr) + lpm(rr)) * ( gradRpm - l(:,ee) ) ...
                - sgnP0u(rr)*abs(z)*F ...
                - sgnP0u(rr)*atanterm(rr) * atanterm(rr) * z(rr) * [0;0;1] );  
    end
    %I1 = I1 + bsxfun(@times,u(:,ee)/2, (R02 .* logterm + lplus.*Rplus - lminus.*Rminus).').'; 
end

%Constructing the final matrix 
%gradVtri = I0*a.' gradI1*a' + bsxfun(@times, b.',I0);