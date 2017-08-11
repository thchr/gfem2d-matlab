function V = calcCoulomb(p,t,method)
%CALL :        V = calcCoulomb(p,t,method)
%Computes the Coulomb matrix associated with a (non-periodic) mesh 
%specified by (p,t) using either of two methods (specified by 'method'):
%   'analytical' : Analytical integration of the linearly varying charge
%                  density weighted by the Coulomb interaction. (default)
%                  More accurate than 'numquad' but takes longer evaluation
%                  See Wilton et al. IEEE Trans. Ant. Prop. 32, 276 (1984)
%   'numquad' : Numerical 16-point quadrature by subdivision of each
%               triangular element into 16 subtriangles, with quadrature
%               performed at their centerpoints.


if ~exist('method','var') || isempty(method)
    method = 'analytical';
end
    
switch method 
    case 'analytical'
        V = calcCoulombAnalytical(p,t);
    case 'numquad'
        V = calcCoulombNumQuad(p,t);
end

