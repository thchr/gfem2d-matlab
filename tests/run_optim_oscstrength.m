clear all; close all; 

%% DEFINE A 'SEED' MESH
type = 'ellipse';
N = 126; 
switch type
    case 'disk'
    theta = linspace(0,2*pi,N+1); theta = theta(1:end-1); 
    periphery(:,1) = cos(theta)/sqrt(pi); 
    periphery(:,2) = sin(theta)/sqrt(pi); 
    case 'ellipse'
    theta = linspace(0,2*pi,N+1); theta = theta(1:end-1); 
    periphery(:,1) = cos(theta)/sqrt(pi)*8; 
    periphery(:,2) = sin(theta)/sqrt(pi);     
end

[os,zeta] = optim_oscstrength(periphery,1,[],1);
%% OPTIMIZATION
%fminunc(@(xy) -max(optim_oscstrength([xy;[flipud(xy(2:N/2,1)),-flipud(xy(2:N/2,2))]],...
%    1,[],1)),periphery(1:N/2+1,:))
initpars = [8.7400,0.3646,-0.8589,2.2944];%[6.4,1,-.64,2]; 
initpars = [38.8435,0.5164,-4.8324,3.2006];
initpars = [53.7888,0.3964,-3.0917,3.6358]
%fminsearch(@(pars) -max(optim_oscstrength(geomEllipseAnnulus(pars(1:2),pars(3:4),N),1,[],1)),initpars)
A = [-1,0,1;
     0,-1,0; 
     0,0,1;
     0,0,0]; 
b = [0.05,0.05,.9,.9].';
%fmincon(@(pars) -max(optim_oscstrength(geomPeanut(pars(1),pars(2),pars(3),N),1,[],1)),initpars,A,b)
fminsearch(@(pars) -max(optim_oscstrength(geomPeanut(pars(1),pars(2),pars(3),pars(4),N),1,[],1)),initpars)