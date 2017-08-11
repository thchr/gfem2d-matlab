function y=genexpint(a,x,expscale)
% generalised exponential integral function [S Bocquet DSTO Oct 2013]
%
%DESCRIPTION COPIED FROM MATHWORKS PAGE:
%This function computes the generalized exponential integral E_a(x)
%for positive real parameter a and argument x. Call it as 
%y=genexpint(a,x) or y=genexpint(a,x,expscale). If the optional third 
%input argument expscale is set to true, the output is exp(x)*E_a(x), 
%which is finite for large x where exp(x) overflows and E_a(x) underflows. 
%The code uses a MATLAB translation of the FORTRAN function DGAMIC from the
%SLATEC library. DGAMIC computes the upper incomplete gamma function for
%negative real parameter, using the algorithm of Gautschi (ACM Trans. 
%Math. Soft. 5(4) pp 466-481, 1979). For x>1, the Legendre continued 
%fraction is used to calculate the function G(1-a)=exp(x)*E_a(x). 
%For x<=1, the generalized exponential integral is obtained from 
%the relationship E_a(x)=x^(a-1)*Gamma(1-a,x). 
%(Note that the incomplete gamma function parameter 1-a can be 
%negative, so the MATLAB function gammainc cannot be used here 
%as it is limited to positive real parameters.)


if any(x(:)<0), error('genexpint:negx','x must be >=0'); end
if nargin<3, expscale=false; end
if isscalar(x), sz=size(a); else sz=size(x); end
if isscalar(x) && ~isscalar(a), x=repmat(x,numel(a),1); a=a(:); end
if isscalar(a) && ~isscalar(x), a=repmat(a,numel(x),1); x=x(:); end
y=NaN(size(x));
b=x>1;
y(b)=glcf1(1-a(b),x(b)); % use Legendre continued fraction for x>1
if ~expscale, y(b)=y(b).*exp(-x(b)); end
s=(x>0 & x<=1);
y(s)=gfn(a(s),x(s),expscale); % use dgamic for 0<x<=1
y(x==0)=Inf;
z1=(x==0 & a>1);
y(z1)=1./(a(z1)-1);
y=reshape(y,sz);
end

function y=gfn(a,x,expscale)
xsc=(a-1).*log(x);
if expscale, xsc=xsc+x; end % scale by exponential
y=exp(xsc).*gamic(1-a,x);
end
