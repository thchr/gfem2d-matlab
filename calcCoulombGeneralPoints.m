function V = calcCoulombGeneralPoints(p,t,r)
%CALL:        V = calcCoulombGeneralPoints(p,t,r)
%DESCRIPTION: Calculates the Coulomb matrix for interaction between a 
%triangular mesh (p,t) and arbitrary (in-plane) points r (which is an Nx2
%array for xy-coords or an Nx3 array for xyz coords).
%Uses the Wilton scheme for maximal accuracy. See reference:
%   Wilton, Rao, Glisson, Schaubert, Al-Bundak, and Butler, 
%   'Potential integrals for uniform and linear source distributions on 
%   polygonal and polyhedral domains', 
%   IEEE Trans. Antennas Propag. (1984), DOI: 10.1109/TAP.1984.1143304.
%and the accompanying note you wrote yourself 'Coulomb integrals over
%triangular 2D mesh elements: applications to FEM-BEM method for graphene
%plasmonics'
%Modified on Dec. 19, 2016 to allow non-in-plane input 'r'.

%If r is 2D (x,y): zeroz = 1. Otherwise, if r is 3D (x,y,z): zeroz = 0.
zeroz = size(r,2) == 2; 

%Number of vertices and elements (triangles)
Nvert = size(p,1);
Ntri = size(t,1);

fprintf('Calculates the Coulomb matrix by analytical evaluation of integral terms (Wilton-scheme) V\n')

%Preallocation
V = zeros(size(r,1),Nvert);

timeV = tic;
for j = 1:Ntri
    if mod(j,500) == 0; fprintf('   j-loop:   %g/%g (%.1f min/%.1f min)\n',j,Ntri,toc(timeV)/60,toc(timeV)/60/(j/Ntri)); end
    l = t(j,1); m = t(j,2); n = t(j,3);
    rl = p(l,:); rm = p(m,:); rn = p(n,:);
    rtri = [rl;rm;rn]; %Vertices of j'th triangle (coordinates as rows [x y; x y; xy])
    
    %Each element is calculated by numerical evaluation of the analytical
    %expressions obtained by Wilton et al., and then added into the total
    %Coulomb matrix at the vertex sites
    if zeroz %"Request points" are xy-only
        V(:,[l m n]) = V(:,[l m n]) + real(calcWiltonContribution(rtri,r)); 
    else     %"Request points" are xy _and_ z
        V(:,[l m n]) = V(:,[l m n]) + real(calcWiltonContribution(rtri,r(:,1:2),r(:,3))); 
    end
    %Small rounding errors in the logarithms (I suspect) may give rise to small
    %imaginary parts (near the machine error): we purge them immediately.
    %(in addition, this keeps V as a real double, which is advantageous for
    %addition - otherwise, the procedure may be significantly slowed)
end


fprintf('\n')

function Vtri = calcWiltonContribution(rtri,r,z)
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
rminus = rtri.';
rplus = circshift(rtri,-1).';
l = bsxfun(@times,(rplus - rminus),1./sqrt(sum((rplus-rminus).^2,1))); %Edge tangential unit vector
u = [l(2,:);-l(1,:)]; %Edge normal unit vector

I0 = zeros(size(r,1),1); I1 = zeros(size(r,1),2); zterm = 0;
for ee = 1:3 %This part of the implementation could probably be vectorized and hence improved; too much effort at the moment
    rplusmr = bsxfun(@minus,rplus(:,ee), r.');
    rminusmr = bsxfun(@minus,rminus(:,ee), r.');
    lplus = rplusmr.'*l(:,ee);
    lminus = rminusmr.'*l(:,ee);
    
    P0u = rplusmr.'*u(:,ee); %A term propto dot(lhat,u) is dropped, as lhat and u are orthog [this is dot(vekP0,u)]
    R02 = P0u.^2;
    if nargin == 2 %"Request points" are xy-only
        Rplus = sqrt( sum( rplusmr.^2, 1 ) ).';
        Rminus= sqrt( sum( rminusmr.^2, 1 ) ).';
    else            %"Request points" are xyz
        R02 = R02 + z.^2;                                 %These are the only modification necessary 
        Rplus = sqrt( sum( rplusmr.^2, 1 )  + z.^2.' ).'; %to the scheme for non-in-plane-coords: re-
        Rminus= sqrt( sum( rminusmr.^2, 1 ) + z.^2.' ).'; %definitions of the capital R variables and
        zterm = abs(z)./abs(P0u) .* ( atan(abs(P0u).*lplus ./(R02 + abs(z).*Rplus )) ... %a new z-dependent  
                                      -atan(abs(P0u).*lminus./(R02 + abs(z).*Rminus)) ); %term in I0
        zterm(z == 0) = 0; 
    end
        
    logterm = log( (Rplus + lplus)./(Rminus + lminus) );
    logterm( isnan(logterm) | isinf(logterm) ) = 0; %Remove singular terms which are later multiplied by zero

    %Adding into the integral terms
    I0 = I0 + P0u .* ( logterm - zterm ); 
    I1 = I1 + bsxfun(@times,u(:,ee)/2, (R02 .* logterm + lplus.*Rplus - lminus.*Rminus).').'; 
end

%Constructing the final matrix 
Vtri = I1*a' + bsxfun(@times, b.',I0);