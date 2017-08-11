clc; clear all; close all;
[cols,nams] = flatcolors;
r1 =  [1,.5];
r2 =  [2.5,1.5];
r3 = [.5,2.5];
rtri = [r1;r2;r3]; 
r =  [4,5];
area = meshArea(rtri,[1,2,3]);

%% Figure
set_figsize([],10,10);
plot(rtri([1:end,1],1),rtri([1:end,1],2),'-','Color',cols{24}); hold on
plot(r(:,1),r(:,2),'o','MarkerFaceColor',cols{1},'MarkerEdgeColor',cols{4})
axis equal
%xlim(minmax([0;rtri(:,1);.3+r(1)]))
%ylim(minmax([0;rtri(:,2);.3+r(2)]))

%%
I0 = 0; 
I1 = [0;0];
for k = 1:3;
   rm = rtri(k,:).';
   rp = circshift(rtri,-1); rp = rp(k,:).';
   l = (rp - rm)/norm(rp-rm);
   rmid = rp; 
   plot(rp(1) + [0,l(1)], rp(2) + [0,l(2)],':','Color',cols{3})
   u = cross([l;0],[0;0;1]); u = u(1:2);
   plot(rp(1) + [0,u(1)], rp(2) + [0,u(2)],':','Color',cols{6})
   lp = (rp - r.').'*l;
   lm = (rm - r.').'*l;
   
   P0 = (rp - r.') - lp*l;
   R0 = norm( (rp -r.').'*u);
   Rp = norm(rp - r.'); Rm = norm(rm - r.');
   
   % I0
   if ~all(abs(P0)<1e-11); 
       I0 = I0 + P0.'*u*log( (Rp + lp)/(Rm + lm) );
   end
   
   %I1
   if abs(R0) > 1e-11
       I1 = I1 + u*( (R0)^2*log( (Rp + lp)/(Rm + lm) ))/2;
   else
       [k, Rp + lp,Rm + lm]
   end
   I1 = I1 + u*(lp*Rp - lm*Rm)/2;
   %I1 = I1 + u*( (R0)^2*log( (Rp + lp)/(Rm + lm) )+ lp*Rp - lm*Rm)/2;
end


rpx = @(eta,xi) (1-eta-xi)*r1(:,1) + eta*r2(:,1) + xi*r3(:,1);
rpy = @(eta,xi) (1-eta-xi)*r1(:,2) + eta*r2(:,2) + xi*r3(:,2);
rp = @(eta,xi) [rpx(eta,xi);rpy(eta,xi)];
invR = @(eta,xi) 1./sqrt( (rpx(eta,xi) - r(1)).^2+(rpy(eta,xi) - r(2)).^2 );

int1 = @(eta) quadgk(@(xi) invR(eta,xi),0,1-eta);
etaspace =0.01; int = 0;
for eta = etaspace/2:etaspace:1-etaspace/2;
    int = int + etaspace*int1(eta);
end
int = 2*area*int;
[I0,int]

invRxpx = @(eta,xi) (rpx(eta,xi) - r(1)).*invR(eta,xi);
invRypy = @(eta,xi) (rpy(eta,xi) - r(2)).*invR(eta,xi);
int1x = @(eta) quadgk(@(xi) invRxpx(eta,xi),0,1-eta);
int1y = @(eta) quadgk(@(xi) invRypy(eta,xi),0,1-eta);
etaspace =0.01; intx = 0; inty = 0;
for eta = etaspace/2:etaspace:1-etaspace/2;
    intx = intx + etaspace*int1x(eta);
    inty = inty + etaspace*int1y(eta);
end
[I1,[intx;inty]*2*area]

%%

rho1 = 1; 
rho2 = 3; 
rho3 = -1;

R = [0,1;-1,0]; 
denom = (R*(r3.' - r1.')).'*(r1-r2).';
peta = -(R*(r3-r1).')/denom;
pxi = -(R*(r1-r2).')/denom;

a(:,1) = -(peta+pxi);   b(1) = 1 - a(:,1).'*(r1-r).';
a(:,2) = peta;          b(2) = -a(:,2).'*(r1-r).';
a(:,3) = pxi;           b(3) = -a(:,3).'*(r1-r).';


V = sum((a.'*I1 + b.'*I0).*[rho1;rho2;rho3]);

rho = @(eta,xi) (1-eta-xi)*rho1 + eta*rho2 + xi*rho3;
int1final = @(eta) quadgk(@(xi) invR(eta,xi).*rho(eta,xi),0,1-eta);
etaspace =0.001; intfinal = 0; 
for eta = etaspace/2:etaspace:1-etaspace/2;
    intfinal = intfinal + etaspace*int1final(eta);
end
[V,intfinal*2*area]