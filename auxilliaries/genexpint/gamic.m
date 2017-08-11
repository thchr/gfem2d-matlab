function y=gamic(a,x)
% upper incomplete gamma function
% uses dgamic for matrix arguments [S Bocquet DSTO Aug 2013]

if isscalar(a)
  y=NaN(size(x));
  for i=1:numel(x)
    y(i)=dgamic(a,x(i));
  end
elseif isscalar(x)
  y=NaN(size(a));
  for i=1:numel(a)
    y(i)=dgamic(a(i),x);
  end
elseif size(a)==size(x)
  y=NaN(size(x));
  for i=1:numel(x)
    y(i)=dgamic(a(i),x(i));
  end
else
  error('gamic:arraysize','incompatible input arrays a, x')
end
end