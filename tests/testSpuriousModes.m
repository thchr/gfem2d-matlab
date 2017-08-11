clear all; close all; clc

addpath(genpath('..')); cols = flatcolors;
%ConstantsUnits0;

N = 75;
geom = 'disk';
switch geom
    case 'disk'
        [p,t] = geomDisk(N,7.5/N,1);
    case 'triangle'
        [p,t] = geomEquiTriangle(1,.0425,1);
    case 'ellipse'
        a = 1.1; b = .9;
        theta = linspace(0,2*pi,N+1).'; theta = theta(1:end-1) + (theta(2)-theta(1))/2;
        nodes = [a*cos(theta),b*sin(theta)];
        [p,t] = geomPolygon(nodes,5/N,1);
    case 'ring'
        rin = .15;
        thetaout = linspace(0,2*pi,N+1).'; thetaout = thetaout(1:end-1) + (thetaout(2)-thetaout(1))/2;
        thetain = linspace(0,2*pi,round(N*sqrt(rin))+1).'; thetain = thetain(1:end-1) + (thetain(2)-thetain(1))/2;
        nodesout = [cos(thetaout),sin(thetaout)]; nodesin = rin*[cos(thetain),sin(thetain)];
        nodes = {nodesout ,nodesin};
        [p,t] = geomPolygon(nodes,5/N,1);
end

%%
[~,DRi,DL] = evalc('calcDifferential(p,t)'); %Run without textual output

[~,DRa] = evalc('calcDifferential(p,t,[0,1;-1,0])'); %Run without textual output
[~,V] = evalc('calcCoulomb(p,t)');

[areas,areavec] = meshArea(p,t);

[eigV0,lambda0]=eig(DL\(DRi*V));
%%
zetac = 0.02;
%[eigV,lambda] = polyeig(1i/(2*pi)*DRa*V,-1/(2*pi)*DRi*V - zetac*DL,0,zetac*DL);
C = [zeros(size(V))             , eye(size(V))                                  , zeros(size(V));  %Companion matrix
    zeros(size(V))             , zeros(size(V))                                , eye(size(V))  ;
    (-1i/(2*pi*zetac))*(DL\(DRa*V)), (1/(2*pi*zetac))*(DL\(DRi*V)) + eye(size(V)), zeros(size(V)) ];
[eigV,lambda]=eig(C,'vector');
%{
%% SVD approach
[Usvd,S,Vsvd] = svd(C); Sn = S; Sn(abs(S)<1e-2) = 0; plot(diag(S),'.','color',cols{14}); hold on; plot(diag(Sn),'o','color',cols{1}); hold off
Csvd = Usvd*Sn*Vsvd'; 
[eigV1,lambda1] = eig(C,'vector');
%}
%% Projector approach
% nareavec = areavec/sqrt(areavec*areavec.');
% P = eye(Nv,Nv) - nareavec.'*nareavec; zeromat = zeros(Nv,Nv);
% C1 =[zeros(size(V))             , eye(size(V))                                  , zeros(size(V));  %Companion matrix
%     zeros(size(V))             , zeros(size(V))                                , eye(size(V))  ;
%     (-1i/(2*pi*zetac))*P*(DL\(DRa*V)), P*(1/(2*pi*zetac))*(DL\(DRi*V)) + P*eye(size(V)), zeros(size(V)) ];
%[eigV1,lambda1] = eig(C1,'vector');
%{
%% Tolerancing approach
C1 = C; C1(abs(C)<1e-7) = 0;
[eigV1,lambda1]=eig(C1,'vector');
%% Constraint guess
zerovec=zeros(size(areavec));onevec=ones(size(areavec));
C1 = [C      , [ areavec, areavec, areavec].' ;
      areavec, areavec, areavec , 0 ];
[eigV1,lambda1]=eig(C1,'vector');
%}
%% Diagonal addition approach / integrated density penalty
     dA = diag(areavec); zV = zeros(size(V));
    penaltyfac = 10/zetac;
    C1 = C - penaltyfac*[dA, zV, zV; zV, zV, zV; zV, zV, zV]; 
    [eigV1,lambda1]=eig(C1,'vector');
%% 
xlims = [0,1]*55.5;  
inds = lambda > xlims(1) & lambda < xlims(2); inds1 = lambda1 > xlims(1) & lambda1 < xlims(2);
try; close(2); end
set_figsize(2,25,25)
plot(real(lambda(inds)),imag(lambda(inds)),'.','color',cols{14}); hold on; 
plot(real(lambda1(inds1)),imag(lambda1(inds1)),'o','color',cols{1}); hold off
xlim(xlims)

%%
intrho=sum(bsxfun(@times,eigV(1:end/3,:),areavec.'),1);
intrho1=sum(bsxfun(@times,eigV1(1:end/3,:),areavec.'),1);

try; close(3); end
set_figsize(3,25,25)
semilogy(real(lambda),abs(intrho),'.','color',cols{14}); hold on; 
semilogy(real(lambda),abs(intrho1),'o','color',cols{1})

