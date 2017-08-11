function V = calcCoulombAnalytical(p,t)
%CALL:        V = calcCoulombAnalytical(p,t)
%DESCRIPTION: Calculates the Coulomb matrix for interaction between a 
%triangular mesh (p,t).
%Uses the Wilton scheme for maximal accuracy. See reference:
%   Wilton, Rao, Glisson, Schaubert, Al-Bundak, and Butler, 
%   'Potential integrals for uniform and linear source distributions on 
%   polygonal and polyhedral domains', 
%   IEEE Trans. Antennas Propag. (1984), DOI: 10.1109/TAP.1984.1143304.
%and the accompanying note you wrote yourself 'Coulomb integrals over
%triangular 2D mesh elements: applications to FEM-BEM method for graphene
%plasmonics'

%Number of vertices and elements (triangles)
Nvert = size(p,1);
Ntri = size(t,1);

fprintf('Calculates the Coulomb matrix by analytical evaluation of integral terms (Wilton-scheme) V\n')

%Preallocation
V = zeros(Nvert,Nvert);

timeV = tic;
for j = 1:Ntri
    if mod(j,500) == 0; fprintf('   j-loop:   %g/%g (%.1f min/%.1f min)\n',j,Ntri,toc(timeV)/60,toc(timeV)/60/(j/Ntri)); end
    l = t(j,1); m = t(j,2); n = t(j,3);
    rl = p(l,:); rm = p(m,:); rn = p(n,:);
    rtri = [rl;rm;rn]; %Vertices of j'th triangle (coordinates as rows [x y; x y; x y])
    
    %Each element is calculated by numerical evaluation of the analytical
    %expressions obtained by Wilton et al., and then added into the total
    %Coulomb matrix at the vertex sites
    V(:,[l m n]) = V(:,[l m n]) + real(calcWiltonContribution(rtri,p)); 
    %Small rounding errors in the logarithms (I suspect) may give rise to small
    %imaginary parts (near the machine error): we purge them immediately.
    %(in addition, this keeps V as a real double, which is advantageous for
    %addition - otherwise, the procedure may be significantly slowed)
end

fprintf('\n')

function Vtri = calcWiltonContribution(rtri,r)
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

I0 = zeros(size(r,1),1); I1 = zeros(size(r,1),2);
for ee = 1:3 %This part of the implementation could probably be vectorized and hence improved; too much effort at the moment
    rplusmr = bsxfun(@minus,rplus(:,ee), r.');
    rminusmr = bsxfun(@minus,rminus(:,ee), r.');
    lplus = rplusmr.'*l(:,ee);
    lminus = rminusmr.'*l(:,ee);
    
    P0u = rplusmr.'*u(:,ee); %A term propto dot(lhat,u) is dropped, as lhat and u are orthog [this is dot(vekP0,u)]
    R02 = P0u.^2;
    Rplus = sqrt( sum( rplusmr.^2, 1 ) ).';
    Rminus =sqrt( sum( rminusmr.^2, 1 ) ).';
    
    logterm = log( (Rplus + lplus)./(Rminus + lminus) );
    logterm( isnan(logterm) | isinf(logterm) ) = 0; %Remove singular terms which are later multiplied by zero
    
    %Adding into the integral terms
    I0 = I0 + P0u .* logterm; 
    I1 = I1 + bsxfun(@times,u(:,ee)/2, (R02 .* logterm + lplus.*Rplus - lminus.*Rminus).').'; 
end

%Constructing the final matrix 
Vtri = I1*a' + bsxfun(@times, b.',I0);