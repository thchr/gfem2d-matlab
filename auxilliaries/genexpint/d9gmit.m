function d9gmitResult=d9gmit(a, x, algap1, sgngam) % ALX deleted as unused
%***BEGIN PROLOGUE  D9GMIT
%***SUBSIDIARY
%***PURPOSE  Compute Tricomi's incomplete Gamma function for small
%            arguments.
%***LIBRARY   SLATEC (FNLIB)
%***CATEGORY  C7E
%***TYPE      DOUBLE PRECISION (R9GMIT-S, D9GMIT-D)
%***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, SMALL X,
%             SPECIAL FUNCTIONS, TRICOMI
%***AUTHOR  Fullerton, W., (LANL)
%***DESCRIPTION
%
% Compute Tricomi's incomplete gamma function for small X.
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
%***end PROLOGUE  D9GMIT

persistent bot ep 

if isempty(bot), bot=log(realmin); end;
if isempty(ep), ep=0.25*eps; end;

if(x <= 0.0d0)
  error('d9gmit:xneg', 'x should be > 0')
end;

ma = fix(a + 0.5d0); % ma is an integer
if(a < 0.0d0)
  ma = fix(a - 0.5d0);
end;
aeps = a - ma;

ae = a;
if(a <(-0.5d0))
  ae = aeps;
end;

t = 1.0d0;
te = ae;
s = t;
for  k=1:200;
  fk = k;
  te = -x.*te./fk;
  t = te./(ae+fk);
  s = s + t;
  if(abs(t) < ep*abs(s)), break; end;
end;
if(abs(t) >= ep*abs(s))
  error('d9gmit:noconv','No convergence in 200 terms of Taylor-s series')
end;

if(a >=(-0.5d0))
  algs = -algap1 + log(s);
  d9gmitResult = exp(algs);
  return;
end;

algs = -dlngam(1.0d0+aeps) + log(s);
s = 1.0d0;
m = -ma - 1;
if(m ~= 0)
  t = 1.0d0;
  for  k=1:m;
    t = x.*t./(aeps-(m+1-k));
    s = s + t;
    if(abs(t) < ep*abs(s)), break; end;
  end;
end;

d9gmitResult = 0.0d0;
algs = -ma.*log(x) + algs;
if(s == 0.0d0 || aeps == 0.0d0)
  d9gmitResult = exp(algs);
  return;
end;

sgng2 = sgngam .* sign(s);
alg2 = -x - algap1 + log(abs(s));

if(alg2 > bot)
  d9gmitResult = sgng2 .* exp(alg2);
end;
if(algs > bot)
  d9gmitResult = d9gmitResult + exp(algs);
end;
return;

end %function d9gmit
