clear all; close all; clc

addpath(genpath('..')); cols = flatcolors;
%ConstantsUnits0;

N = 30;
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
[DRi,DL] = calcDifferentialProjection(p,t); %Run without textual output
[DRio] = calcDifferential(p,t);

try; close(5); end; set_figsize(5,25,25);
imagesc(full(DRi)); colormap(bluewhitered_mod); colorbar; drawnow

[~,V] = evalc('calcCoulomb(p,t)');
%%
areas = meshArea(p,t); Nv = size(p,1);
areavec = zeros(1,Nv);
for jj = 1:size(t,1)
    areavec(t(jj,:)) = areavec(t(jj,:)) + areas(jj)/3;
end
nareavec = areavec/sqrt(areavec*areavec.');
P = eye(Nv,Nv) - nareavec.'*nareavec; zeromat = zeros(Nv,Nv);

[eigV,zeta]=eig(DL\(DRi*V)/(2*pi),'vector');
[eigV1,zeta1]=eig(DL\(DRio*V)/(2*pi),'vector');

%%

xlims = [-1e-3,1]*10;  
inds = zeta > xlims(1) & zeta < xlims(2); inds1 = zeta1 > xlims(1) & zeta1 < xlims(2);
try; close(2); end
set_figsize(2,25,25)
plot(real(zeta(inds)),imag(zeta(inds)),'.','color',cols{14}); hold on; 
plot(real(zeta1(inds1)),imag(zeta1(inds1)),'o','color',cols{1}); hold off
xlim(xlims)
%%
intrho=sum(bsxfun(@times,eigV,areavec.'),1);
intrho1=sum(bsxfun(@times,eigV1,areavec.'),1);

try; close(3); end
set_figsize(3,15,15)
semilogy(real(zeta),abs(intrho),'.','color',cols{14}); hold on; 
semilogy(real(zeta1),abs(intrho1),'o','color',cols{1}); hold off; 