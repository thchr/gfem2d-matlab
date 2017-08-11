clear all; close all; clc

addpath(genpath('..'))
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

DRaeff = DL\DRa;
%%
%figure
%surf(log(abs(DL\DRa))); colormap(hot); colorbar; view(2); %zlim([-10,10]); caxis([-10,10])
%%
[phi,r]=cart2pol(p(:,1),p(:,2));
f = r;%.*exp(1i*phi);
Df = DL\(DRa*f);
set_figsize([],24,20);
trisurf(t,p(:,1),p(:,2),log10(abs(Df))); view(2); colormap(hot); colorbar; %axis equal