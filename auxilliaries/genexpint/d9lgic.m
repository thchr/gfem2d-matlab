function d9lgicResult=d9lgic(a, x, alx)
%***BEGIN PROLOGUE  D9LGIC
%***SUBSIDIARY
%***PURPOSE  Compute the log complementary incomplete Gamma function
%            for large X and for A .LE. X.
%***LIBRARY   SLATEC (FNLIB)
%***CATEGORY  C7E
%***TYPE      DOUBLE PRECISION (R9LGIC-S, D9LGIC-D)
%***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, LARGE X,
%             LOGARITHM, SPECIAL FUNCTIONS
%***AUTHOR  Fullerton, W., (LANL)
%***DESCRIPTION
%
% Compute the log complementary incomplete gamma function for large X
% and for A .LE. X.
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
%***END PROLOGUE  D9LGIC

persistent ep 

if isempty(ep),   ep =0.25.*eps;  end;

xpa = x + 1.0d0 - a;
xma = x - 1.0d0 - a;

r = 0.0d0;
p = 1.0d0;
s = p;
for  k=1:300;
  fk = k;
  t = fk.*(a-fk).*(1.0d0+r);
  r = -t./((xma+2.0d0.*fk).*(xpa+2.0d0.*fk)+t);
  p = r.*p;
  s = s + p;
  if(abs(p) < ep*s), break; end;
end;
if(abs(p) >= ep*s)
  error('d9lgic:noconv','No convergence in 300 terms of continued fraction')
end;

d9lgicResult = a.*alx - x + log(s./xpa);

return;
end %function d9lgic
