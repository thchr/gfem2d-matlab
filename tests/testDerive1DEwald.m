clear all; clc; close all;
cols=flatcolors();

R = [1,0]; 
G = [1,0]*2*pi;
k = .2*G;
kR = k*R(:);

Nr = 100;
rp = [0,0]; r(:,1) = linspace(-1,1,Nr); r(:,2) = linspace(.3,.3,Nr);
dr=bsxfun(@minus,r,rp);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmaxd = 1500000; 
stepn = 15000;
nlist = nmaxd:-stepn:1; if nlist(end)~= 1; nlist(end+1) = 1; end
tic
Wd = zeros(size(r,1),size(rp,1));
for jj = 1:numel(nlist)-1; %In chunks to avoid overloading available memory
    list = nlist(jj)-1:-1:nlist(jj+1);
    %negative n
    Wd = Wd  +  sum(bsxfun(@rdivide,exp(-1i*kR*list),pdist2(r,bsxfun(@plus,rp,bsxfun(@times,-list.',R)))),2); 
    %positive n
    Wd = Wd  +  sum(bsxfun(@rdivide,exp(1i*kR*list),pdist2(r,bsxfun(@plus,rp,bsxfun(@times,list.',R)))),2);
    if mod(jj,25)==0;fprintf('%g/%g\n',jj,numel(nlist)-1); end
end
Wd = Wd + 1./pdist2(r,rp); %n=0 term
toc
  

Wd = Wd.*exp(-1i*(k(1)*dr(:,1)+k(2)*dr(:,2)));
Wd=Wd.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
E = sqrt(pi)/norm(R,2);
nmax = 4; 

W1 = zeros(size(r,1),1);
% for n=-nmax:nmax
%     kGn = k-G*n;
%     abskGn = norm(kGn,2);
%     W1 = W1 + exp(1i*(kGn(1)*dr(:,1) + kGn(2)*dr(:,2)))./abskGn.*erfc(abskGn/(2*E));
% end
% W1 = 2*pi/norm(R,2)*exp(-1i*(k(1)*dr(:,1)+k(2)*dr(:,2))).*W1;
for n=-nmax:nmax
     kGn = k+G*n;
     abskGn = norm(kGn,2);
     %intterm = quadgk(@(x) 2*cos( (k(1)+G(1)*n)*x ).*erf(x*E)./x,0,1000,'AbsTol',1e-10,'RelTol',1e-10,'MaxInterValCount',10000);
     %intterm = quadgk(@(x) 2*cos( (k(1)+G(1)*n)*x ).*erf(sqrt(x.^2+dr(1,2).^2)*E)./sqrt(x.^2+dr(1,2).^2),0,1000,'AbsTol',1e-10,'RelTol',1e-10,'MaxInterValCount',10000);
     intterm = quadgk(@(s) 1./s.*exp(-(k(1)+G(1)*n).^2/(4*E^2)./s - E^2*dr(1,2)^2*s),0,1,'AbsTol',1e-10,'RelTol',1e-10,'MaxInterValCount',10000);
     %intterm = 1/norm(R,2)*expint(1/4* (k(1)+G(1)*n)^2/E^2);
     W1 = W1 + exp(1i*(G(1)*n*dr(:,1) + G(2)*n*dr(:,2))).*intterm;
end


W2 = zeros(size(r,1),1);
for n=-nmax:nmax
    absrRn = pdist2(dr,n*R); 
    W2 = W2 + exp(1i*kR*n)./absrRn.*erfc(absrRn*E);
end
W2 = exp(-1i*(k(1)*dr(:,1)+k(2)*dr(:,2))).*W2;


We = W1+W2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
figure
semilogy(r(:,1),abs(Wd),'-','color',cols{1})
hold on
semilogy(r(:,1),abs(We),'--','color',cols{2})
semilogy(r(:,1),abs(W1),':','color',cols{3})
semilogy(r(:,1),abs(W2),':','color',cols{4})
hold off

figure 
plot(r(:,1),Wd-W2.'); 
hold on
plot(r(:,1),abs(real(W1)),':','color',cols{3})
hold off

figure 
semilogy(r(:,1),Wd-We.'); 

