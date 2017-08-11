function d9lgitResult=d9lgit(a, x, algap1)
%***BEGIN PROLOGUE  D9LGIT
%***SUBSIDIARY
%***PURPOSE  Compute the logarithm of Tricomi's incomplete Gamma
%            function with Perron's continued fraction for large X and
%            A .GE. X.
%***LIBRARY   SLATEC (FNLIB)
%***CATEGORY  C7E
%***TYPE      DOUBLE PRECISION (R9LGIT-S, D9LGIT-D)
%***KEYWORDS  FNLIB, INCOMPLETE GAMMA FUNCTION, LOGARITHM,
%             PERRON'S CONTINUED FRACTION, SPECIAL FUNCTIONS, TRICOMI
%***AUTHOR  Fullerton, W., (LANL)
%***DESCRIPTION
%
% Compute the log of Tricomi's incomplete gamma function with Perron's
% continued fraction for large X and for A .GE. X.
%
%***REFERENCES  (NONE)
%***ROUTINES CALLED  D1MACH, XERMSG
%***REVISION HISTORY  (YYMMDD)
%   770701  DATE WRITTEN
%   890531  Changed all specific intrinsics to generic.  (WRB)
%   890531  REVISION DATE from Version 3.2
%   891214  Prologue converted to Version 4.0 format.  (BAB)
%   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
%   900720  Routine changed from user-callable to subsidiary.  (WRB)
% Matlab translation [S Bocquet DSTO Aug 2013]
%***END PROLOGUE  D9LGIT

persistent ep sqeps 

if isempty(ep), ep=0.25*eps; end;
if isempty(sqeps), sqeps=sqrt(eps); end;

if(x <= 0.0d0 || a < x)
  error('d9lgit:xvalue','x should be > 0.0 and <= a')
end;

ax = a + x;
a1x = ax + 1.0d0;
r = 0.0d0;
p = 1.0d0;
s = p;
for  k=1:200;
  fk = k;
  t =(a+fk).*x.*(1.0d0+r);
  r = t./((ax+fk).*(a1x+fk)-t);
  p = r.*p;
  s = s + p;
  if(abs(p) < ep*s), break; end;
end;
if(abs(p) >= ep*s)
  error('d9lgit:noconv','No convergence in 200 terms of continued fraction')
end;

hstar = 1.0d0 - x.*s./a1x;
if(hstar < sqeps)
  warning('d9lgit:lowprec','Result less than half precision')
end;

d9lgitResult = -x - algap1 - log(hstar);
return;

end %function d9lgit
