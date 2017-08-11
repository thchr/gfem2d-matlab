clear all;close all;clc;
cols=flatcolors;
R{1} = [1;.99998]; R{2} = [0;1];
[G,A_uc] = calcReciprocal(R);
dr{1} = .11; dr{2} = .1;
k = [pi&2,pi/3]/inf;

nmax = 1:10;
for ii = 1:numel(nmax);
[WE(ii)] = calcLatticeEwaldCoulomb(dr,k,R,G,nmax(ii));
if ii > 2
semilogy(nmax(1:ii-1),abs((WE(1:ii-1)-WE(ii))./WE(ii)),'-s','color',cols{1})
end
drawnow;
end

%%
figure
loglog(nmax(1:end-1),abs((WE(1:end-1)-WE(end))./WE(end)),'-s','color',cols{1})


%{
nmax = 1:10;
for ii = 1:numel(nmax);
[WE(ii),W(ii),dW(ii)] = calcLatticeEwaldCoulomb(dr,k,R,G,nmax(ii));
subaxis(1,3,1)
plot(nmax(1:ii),real(WE),'-s','color',cols{1});hold on;
plot(nmax(1:ii),real(W),'-s','color',cols{2});hold off
subaxis(1,3,2)
plot(nmax(1:ii),imag(WE),':o','color',cols{1});hold on
plot(nmax(1:ii),imag(W),':o','color',cols{2});hold off;
drawnow;
subaxis(1,3,3)
plot(nmax(1:ii),abs(dW./WE),':o','color',cols{1});hold on
drawnow;
end
%}