function W = calcLatticeDirectSumCoulomb(dr,k,R,nmax)
%Calculate the lattice sum of the effective Coulomb potential in a periodic
%structure characterized by 2D lattice vectors R{1} and R{2} for a nonzero
%momentum k for seperations dr = r - r', for points r and r' within the
%unit cell. The following variable form is expected:
%       dr : 2-element cell array of vector/matrix arrays
%       k  : [1 x 2] row vector
%       R  : 2-element cell array with entries R{i} : [1 x 2] row vector
%       Nmax : positive integer (cutoff in lattice summation)
%EXAMPLE
%     x = linspace(-.5,.5,50); y = x; [X,Y] = meshgrid(x,y);
%     x0 = 0; y0 = 0;
%     dr = { X-x0, Y-y0 };
%     R = { [1,0], [0,1] };
%     k = 2*pi*[0,.25];
%     W = calcLatticeCoulomb(dr,k,R,75)
%     contourf(abs(W)); axis equal; colormap(hot)

%Default Nmax
if ~exist('nmax','var') || isempty(nmax)
    nmax = 75;
end

W = zeros(size(dr{1})); %Preallocate
for n = -nmax:nmax %Sum over lattice vectors
    for m = -nmax:nmax
        W = W + exp(1i*k*(n*R{1}+m*R{2}).')./sqrt( (dr{1}-n*R{1}(1)-m*R{2}(1)).^2 + (dr{2}-n*R{1}(2)-m*R{2}(2)).^2 );        
    end
end

W = exp(-1i * ( k(1)*dr{1} + k(2)*dr{2} )) .* W; %Multiply phase factor
