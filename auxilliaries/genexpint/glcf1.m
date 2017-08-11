function y=glcf1(a,x)
% calculates the function G(a,x)=exp(x)*x^(-a)*Gamma(a,x) using the Legendre
% continued fraction. Based on D9LGIC from SLATEC library.
% [S Bocquet DSTO Oct 2013]
if isscalar(a)
  y=NaN(size(x));
  for i=1:numel(x)
    y(i)=lcf(a,x(i));
  end
elseif isscalar(x)
  y=NaN(size(a));
  for i=1:numel(a)
    y(i)=lcf(a(i),x);
  end
elseif size(a)==size(x)
  y=NaN(size(x));
  for i=1:numel(x)
    y(i)=lcf(a(i),x(i));
  end
else
  error('gamic:arraysize','incompatible input arrays a, x')
end
end

function y=lcf(a,x)

persistent ep 
if isempty(ep),   ep =0.25.*eps;  end

xpa = x + 1.0 - a;
xma = x - 1.0 - a;

r = 0.0;
p = 1.0;
s = p;
for  fk=1:300
  t = fk.*(a-fk).*(1.0+r);
  r = -t./((xma+2.0.*fk).*(xpa+2.0.*fk)+t);
  p = r.*p;
  s = s + p;
  if abs(p) < ep*s, break; end
end
if abs(p) >= ep*s
  error('lcf:noconv','No convergence in 300 terms of continued fraction')
end
y = s./xpa;
end
