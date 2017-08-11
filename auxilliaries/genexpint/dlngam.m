function [dlgam, sgngam]=dlngam(x)
% log(abs(gamma(x)))
% Replaces DLGAMS and DLNGAM in dgamic [S Bocquet DSTO Aug 2013]
dlgam=zeros(size(x)); sgngam=ones(size(x));
p=x>=0;
dlgam(p)=gammaln(x(p));
sinpix=sin(pi*x(~p));
dlgam(~p)=log(pi)-log(abs(sinpix))-gammaln(1-x(~p));
sgngam(~p)=sign(sinpix); % sign of gamma function
ni=(~p & x==floor(x)); % correction for x a negative integer
dlgam(ni)=Inf;
sgngam(ni)=1; % actually sign is indeterminate, but convenient to set to 1
end