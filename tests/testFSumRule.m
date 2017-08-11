clear all; close all; clc;

addpath(genpath('..'))
ConstantsUnits0;

%% SETUP
Neigs = 150;
N = 150;
Nab = 3;
maxab = 1;
alist = 1./sqrt(linspace(pi,1/maxab^2,Nab));%linspace(sqrt(1/pi),maxab,Nab);
alist = 2*[sqrt(linspace(1/(pi*maxab)^2,1/pi,Nab)) alist(2:end)];%[1./(pi*linspace(maxab,sqrt(1/pi),Nab)) alist(2:end)];

%% CALCULATION
ellipse = cell(numel(alist));
for aa = 1:numel(alist);
    a = alist(aa); b = 1/(pi*a); blist(aa) = b; elist(aa) = sign(a-b)*sqrt(1-(min(a,b)./max(a,b)).^2);
    
    theta = linspace(0,2*pi,N+1).'; theta = theta(1:end-1) + (theta(2)-theta(1))/2;
    nodes = [a*cos(theta),b*sin(theta)];
    [~,mesh.p,mesh.t] = evalc('geomPolygon(nodes,8/N,1)');
    
    fprintf('GEOMETRY %g/%g: b/a = %.2f\n',aa,numel(alist),b/a)
    [alpha,eigstruct] = calcPolarizability(mesh,Neigs);%round(size(mesh.p,1)/2));
    
    %Store the data in a cell array of structures
    ellipse{aa}.a = a; ellipse{aa}.b = b;
    ellipse{aa}.bndry = [nodes; nodes(1,:)];
    ellipse{aa}.mesh = mesh;
    ellipse{aa}.alpha = alpha;
    ellipse{aa}.eigstruct = eigstruct;
    ellipse{aa}.N = N;

    fprintf('\n   SUM RULE CHECK (%g lowest modes): x = %.4f, y = %.4f\n\n',Neigs,sum(eigstruct.oscstrength(:,1)),sum(eigstruct.oscstrength(:,2)))
end

