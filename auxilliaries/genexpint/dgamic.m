function dgamicResult=dgamic(a, x)
%***BEGIN PROLOGUE  DGAMIC
%***PURPOSE  Calculate the complementary incomplete Gamma function.
%***LIBRARY   SLATEC (FNLIB)
%***CATEGORY  C7E
%***TYPE      DOUBLE PRECISION (GAMIC-S, DGAMIC-D)
%***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB,
%             SPECIAL FUNCTIONS
%***AUTHOR  Fullerton, W., (LANL)
%***DESCRIPTION
%
%   Evaluate the complementary incomplete Gamma function
%
%   DGAMIC = integral from X to infinity of EXP(-T) * T**(A-1.)  .
%
%   DGAMIC is evaluated for arbitrary real values of A and for non-
%   negative values of X (even though DGAMIC is defined for X .LT.
%   0.0), except that for X = 0 and A .LE. 0.0, DGAMIC is undefined.
%
%   DGAMIC, A, and X are DOUBLE PRECISION.
%
%   A slight deterioration of 2 or 3 digits accuracy will occur when
%   DGAMIC is very large or very small in absolute value, because log-
%   arithmic variables are used.  Also, if the parameter A is very close
%   to a negative INTEGER (but not a negative integer), there is a loss
%   of accuracy, which is reported if the result is less than half
%   machine precision.
%
%***REFERENCES  W. Gautschi, A computational procedure for incomplete
%                 gamma functions, ACM Transactions on Mathematical
%                 Software 5, 4 (December 1979), pp. 466-481.
%               W. Gautschi, Incomplete gamma functions, Algorithm 542,
%                 ACM Transactions on Mathematical Software 5, 4
%                 (December 1979), pp. 482-489.
%***ROUTINES CALLED  D1MACH, D9GMIC, D9GMIT, D9LGIC, D9LGIT, DLGAMS,
%                    DLNGAM, XERCLR, XERMSG
%***REVISION HISTORY  (YYMMDD)
%   770701  DATE WRITTEN
%   890531  Changed all specific intrinsics to generic.  (WRB)
%   890531  REVISION DATE from Version 3.2
%   891214  Prologue converted to Version 4.0 format.  (BAB)
%   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
%   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
% Matlab translation [S Bocquet DSTO Aug 2013]
%***end PROLOGUE  DGAMIC

persistent alneps bot ep sqeps 

if isempty(alneps), alneps=-log(0.5*eps); end;
if isempty(bot), bot=log(realmin); end;
if isempty(ep), ep=0.25*eps; end;
if isempty(sqeps), sqeps=sqrt(eps); end;

if(x < 0.0d0)
  error('dgamic:xneg','x is negative')
elseif(x == 0.0d0)
  if(a <= 0.0d0)
    %     error('dgamic:undef','x = 0 and a <= 0 so dgamic is undefined')
    dgamicResult = NaN;
  else
    dgamicResult = exp(dlngam(a+1.0d0) - log(a));
  end
  return
end
alx = log(x);
sga = 1.0d0;
if(a ~= 0.0d0)
  sga = sign(a);
end
ainta = fix(a + 0.5d0.*sga);
aeps = a - ainta;

izero = false;
if(x < 1.0d0)
  if (a <= 0.5d0 && abs(aeps) <= 0.001d0)
    e = 2.0d0;
    if((-ainta) > 1.0d0)
      e = 2.0d0.*(-ainta+2.0d0)./(ainta.*ainta-1.0d0);
    end
    e = e - alx.*x.^(-0.001d0);
    if(e*abs(aeps) <= ep)
      dgamicResult = d9gmic(a, x, alx);
      return
    end
  end
  [algap1, sgngam]=dlngam(a+1.0d0);
  gstar = d9gmit(a, x, algap1, sgngam);  % alx deleted as unused
  if(gstar == 0.0d0)
    izero = true;
  else
    alngs = log(abs(gstar));
    sgngs = sign(gstar);
  end
else
  if(a < x)
    dgamicResult = exp(d9lgic(a, x, alx));
    return
  end
  sgngam = 1.0d0;
  algap1 = dlngam(a+1.0d0);
  sgngs = 1.0d0;
  alngs = d9lgit(a, x, algap1);
end
% EVALUATION OF DGAMIC(A,X) in TERMS OF TRICOMI-S INCOMPLETE GAMMA FN.
h = 1.0d0;
if(~izero)
  t = a.*alx + alngs;
  if(t > alneps)
    sgng = -sgngs.*sga.*sgngam;
    t = t + algap1 - log(abs(a));
%     if(t < bot)
%       xerclr;
%     end;
    dgamicResult = sgng.*exp(t);
    return
  end
  if(t >(-alneps))
    h = 1.0d0 - sgngs.*exp(t);
  end
%   if(abs(h) < sqeps)
%     xerclr;
%   end;
  if(abs(h) < sqeps)
    warning('dgamic:lowprec','Result less than half precision')
  end
end
sgng = sign(h).*sga.*sgngam;
t = log(abs(h)) + algap1 - log(abs(a));
% if(t < bot)
%   xerclr;
% end;
dgamicResult = sgng.*exp(t);
return

end %function dgamic
