clc; clear all; close all; 
cols=flatcolors;
addpath(genpath('..'))

%Inclusion boundary
adivd = 2;
Ncirc = 12; 
theta = linspace(0,2*pi,Ncirc+1); theta = theta(1:end-1);
inclusion = [];1/(2*adivd)*[cos(theta);sin(theta)].';

blochmesh=geomRibbonInclusion('triangular',[1,1],inclusion,25,[],'remesh');


k = [0.35,0]*2*pi;

%Unwrapping the necessary elements of the blochmesh structure
areas = blochmesh.area;
R = blochmesh.Rrib; %Lattice vector
G = blochmesh.Grib; %Reciprocal lattice vector
remesh = blochmesh.remesh;

%If unspecified, choose default value
if ~exist('nmax','var') || isempty(nmax);
    nmax = 5 + (norm(k) == 0)*2;
end
loop = -nmax:nmax; 

nn=20;
r = blochmesh.remesh.p; rp = [-.6,-.3];blochmesh.remesh.p(nn,:);

%Unit cell area
L_uc = norm(R,2);

%Direct and reciprocal lattice vectors  needed in the Ewald scheme are
%computed once and for all here
Rn = bsxfun(@times,loop',R); %[(2*nmax+1) x 2] vector
Gn = bsxfun(@times,loop',G); %[(2*nmax+1) x 2] vector
tic
%The optimal integral splitting parameter in the Ewald scheme (guessing)
E = sqrt(pi)/L_uc;

%Reocurring exponential elements needed in the Ewald scheme
expiGnr = exp(-1i*( bsxfun(@times, r(:,1), Gn(:,1).') + bsxfun(@times, r(:,2), Gn(:,2).') ) );
expiRnk = exp(1i*(k(1)*Rn(:,1) + k(2)*Rn(:,2)));
expikr =   exp(-1i*( k(1)*r(:,1)+k(2)*r(:,2) ) );

%Reoccuring erfc elements needed in Ewald summation scheme
abskGn = sqrt( (k(1) - Gn(:,1)).^2 + (k(2) - Gn(:,2)).^2 );
erfcabskGn = erfc(abskGn/(2*E))./abskGn;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%     SUBFUNTIONS     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%The rp-dependent part of the nm-sum in matrix form (vector operation)
expiGnrp = exp(1i*( bsxfun(@times, rp(:,1), Gn(:,1).') + bsxfun(@times, rp(:,2), Gn(:,2).' ) ) );
erfcexpGnrp = bsxfun(@times, expiGnrp, erfcabskGn.');

%Calculate W1 by a vectorized approach (faster than looping)
excludekzero = Gn(:,1) == k(1) & Gn(:,2) == k(2); 
if any(excludekzero) %Skip any nm components where Gnm = k , which gives the divergent part - importantly, is dr-independent and so may be excluded
    erfcexpGnrp(:,excludekzero) = []; %Just remove the "offending" column
    expiGnr(:,excludekzero) = [];
end
W1 = erfcexpGnrp*expiGnr.'; %Efficient summation approach by matrix multiplication

W1 = W1*(2*pi/L_uc); %No reshaping necessary




%Number of points in r and rp (assumed in column form)
numr = size(r,1); numrp = size(rp,1);
%Calculate W2 by vector operations only (the following two calculations, i.e.
%absdrRnm and W, take >95% of all evaluation-time for calcLatticeCoulombOptimized!)
absdrRn = pdist2([reshape(bsxfun(@minus, r(:,1).', rp(:,1)),numrp*numr,1),... %Calculates all distances between (r, rp), and Rnm 
                   reshape(bsxfun(@minus, r(:,2).', rp(:,2)),numrp*numr,1)],...%This implementation is faster than a sqrt-implementation
                   Rn);                                                       %because it utilizes a .mex function (so it is C I guess)
W2 = (erfc(absdrRn*E)./absdrRn) * expiRnk; %The matrix multiplication takes care of a summation-type operation
prefactor = bsxfun(@times, expikr.', exp(1i*(k(1)*rp(:,1) + k(2)*rp(:,2))));
W2 =  prefactor(:).*W2;

%Reshape to the size we would get from bsxfun'ing r.' with rp
W2 = reshape(W2, [numrp, numr]);

%%%%%%%%%%% SUMMING TERMS %%%%%%%%%%
W = W1+W2;
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
nmaxd = 5;15000000;1000000; 

kR = k*R(:);
stepn = 1;1500;
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
  
dr=bsxfun(@minus,r,rp);
Wd = Wd.*exp(-1i*(k(1)*dr(:,1)+k(2)*dr(:,2)));
Wd=Wd.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

figure
%trisurf(blochmesh.remesh.t,r(:,1),r(:,2),(real(W)-real(Wd))./real(W),'EdgeColor',cols{1},'FaceColor','none');
trisurf(blochmesh.remesh.t,r(:,1),r(:,2),(real(W)-0*min(min(abs(W)))),'EdgeColor',cols{1},'FaceColor','none');
hold on
trisurf(blochmesh.remesh.t,r(:,1),r(:,2),(real(Wd)-0*min(min(abs(Wd)))),'EdgeColor',cols{2},'FaceColor','none');
hold off
