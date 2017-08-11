function V = calcCoulombNumQuad(p,t)

%Number of vertices and elements (triangles)
Nvert = size(p,1);
Ntri = size(t,1);

%Triangular element areas
areas = meshArea(p,t);

%Gaussian quadrature weights
GaussQuadMat = [ 10, 7, 8, 7, 4, 5, 4, 5, 4,  1, 2, 1, 2, 1, 2,  1 ; ...
                  1, 4, 2, 1, 7, 5, 4, 2, 1, 10, 8, 7, 5, 4, 2,  1 ; ...
                  1, 1, 2, 4, 1, 2, 4, 5, 7,  1, 2, 4, 5, 7, 8, 10 ];
              
fprintf('Calculates the Coulomb matrix via 16-point numerical quadrature V\n')

%Preallocation
V = zeros(Nvert,Nvert);

timeV = tic;
for j = 1:Ntri
    if mod(j,500) == 0; fprintf('   j-loop:   %g/%g (%.1f min/%.1f min)\n',j,Ntri,toc(timeV)/60,toc(timeV)/60/(j/Ntri)); end
    s = areas(j);
    l = t(j,1); m = t(j,2); n = t(j,3);
    rl = p(l,:); rm = p(m,:); rn = p(n,:);
    
    cpnt = calcQuadraturePoints(rl,rm,rn); %Center points for 16 subtriangles in the j'th triangle
       
    % Calculate the inverse distances between all vertex points, p, 
    % and the jth triangle's subtriangle centerpoints cpnt
    invdcn = 1./sqrt( bsxfun(@minus, p(:,1).', cpnt(:,1)).^2 + bsxfun(@minus, p(:,2).', cpnt(:,2)).^2 );

    % Each element is finally the summation of contributions from all sub-
    % triangles, added into the total Coulomb matrix at the vertex sites   
    V(:,[l m n]) = V(:,[l m n]) + (GaussQuadMat*invdcn).'*s;
    
end

V=V/192; %Normalize by (straggling numerical factors)
fprintf('\n')