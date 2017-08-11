function d9gmicResult=d9gmic(a, x, alx)
%***BEGIN PROLOGUE  D9GMIC
%***SUBSIDIARY
%***PURPOSE  Compute the complementary incomplete Gamma function for A
%            near a negative integer and X small.
%***LIBRARY   SLATEC (FNLIB)
%***CATEGORY  C7E
%***TYPE      DOUBLE PRECISION (R9GMIC-S, D9GMIC-D)
%***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, SMALL X,
%             SPECIAL FUNCTIONS
%***AUTHOR  Fullerton, W., (LANL)
%***DESCRIPTION
%
% Compute the complementary incomplete gamma function for A near
% a negative integer and for small X.
%
%***REFERENCES  (NONE)
%***ROUTINES CALLED  D1MACH, DLNGAM, XERMSG
%***REVISION HISTORY  (YYMMDD)
%   770701  DATE WRITTEN
%   890531  Changed all specific intrinsics to generic.  (WRB)
%   890911  Removed unnecessary intrinsics.  (WRB)
%   890911  REVISION DATE from Version 3.2
%   891214  Prologue converted to Version 4.0 format.  (BAB)
%   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
%   900720  Routine changed from user-callable to subsidiary.  (WRB)
% Matlab translation [S Bocquet DSTO Aug 2013]
%***end PROLOGUE  D9GMIC

persistent bot ep euler 

if isempty(bot), bot=log(realmin); end;
if isempty(ep), ep=0.25*eps; end;
if isempty(euler), euler =0.57721566490153286060651209008240d0;  end;

if(a > 0.0d0)
  error('d9gmic:avalue','a must be near a negative integer')
end;
if(x <= 0.0d0)
  error('d9gmic:xneg', 'x must be > 0')
end;

m = fix(-(a - 0.5d0)); % m is an integer
fm = m;

te = 1.0d0;
t = 1.0d0;
s = t;
for  k=1:200;
  fkp1 = k + 1;
  te = -x.*te./(fm+fkp1);
  t = te./fkp1;
  s = s + t;
  if(abs(t) < ep*s), break; end;
end;
if(abs(t) >= ep*s)
  error('d9gmic:noconv','No convergence in 200 terms of continued fraction')
end;

d9gmicResult = -alx - euler + x.*s./(fm+1.0d0);
if(m == 0), return; end

if(m == 1)
  d9gmicResult = -d9gmicResult - 1.0d0 + 1.0d0./x;
  return
end

te = fm;
t = 1.0d0;
s = t;
mm1 = m - 1;
for  k=1:mm1;
  fk = k;
  te = -x.*te./fk;
  t = te./(fm-fk);
  s = s + t;
  if(abs(t) < ep*abs(s)), break; end;
end;

for  k=1:m;
  d9gmicResult = d9gmicResult + 1.0d0./k;
end;

sgng = 1.0d0;
if(rem(m,2) == 1)
  sgng = -1.0d0;
end;
alng = log(d9gmicResult) - dlngam(fm+1.0d0);

d9gmicResult = 0.0d0;
if(alng > bot)
  d9gmicResult = sgng .* exp(alng);
end;
if(s ~= 0.0d0)
  d9gmicResult = d9gmicResult + (abs(exp(-fm.*alx+log(abs(s)./fm))).*sign(s));
end;

if(d9gmicResult == 0.0d0 && s == 0.0d0)
  warning('d9gmic:underflow','Result underflows')
end;
return;

end %function d9gmic
