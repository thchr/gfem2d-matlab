clear all; close all; clc;

addpath(genpath('..'))
ConstantsUnits0;

%% SETUP
N = 75;
Nab = 25;
maxab = 2;
alist = 1./sqrt(linspace(pi,1/maxab^2,Nab));%linspace(sqrt(1/pi),maxab,Nab);
alist = [sqrt(linspace(1/(pi*maxab)^2,1/pi,Nab)) alist(2:end)];%[1./(pi*linspace(maxab,sqrt(1/pi),Nab)) alist(2:end)];

%% CALCULATION
ellipse = cell(numel(alist));
for aa = 1:numel(alist);
    a = alist(aa); b = 1/(pi*a); blist(aa) = b; elist(aa) = sign(a-b)*sqrt(1-(min(a,b)./max(a,b)).^2);
    
    theta = linspace(0,2*pi,N+1).'; theta = theta(1:end-1) + (theta(2)-theta(1))/2;
    nodes = [a*cos(theta),b*sin(theta)];
    [~,mesh.p,mesh.t] = evalc('geomPolygon(nodes,5/N,1)');
    
    [alpha,eigstruct] = calcPolarizability(mesh);
    
    %Store the data in a cell array of structures
    ellipse{aa}.a = a; ellipse{aa}.b = b;
    ellipse{aa}.bndry = [nodes; nodes(1,:)];
    ellipse{aa}.mesh = mesh;
    ellipse{aa}.alpha = alpha;
    ellipse{aa}.eigstruct = eigstruct;
    ellipse{aa}.N = N;
end

%% SORTING

%Sorting of zeta and oscstrength for the dipole modes along x and y
for aa = 1:numel(alist)
    
    
    if alist(aa) == blist(aa) %The disk needs special care due to symmetry (no natural axis along x,y)
        maxX.oscstrength(aa) = sum(ellipse{aa}.eigstruct.oscstrength(1,:));
        maxY.oscstrength(aa) = sum(ellipse{aa}.eigstruct.oscstrength(2,:));
        Ix = 1; Iy = 2;
    else
        [maxX.oscstrength(aa),Ix] = max(ellipse{aa}.eigstruct.oscstrength(:,1));
        [maxY.oscstrength(aa),Iy] = max(ellipse{aa}.eigstruct.oscstrength(:,2));
    end
    maxX.zeta(aa) = ellipse{aa}.eigstruct.zeta(Ix);
    maxY.zeta(aa) = ellipse{aa}.eigstruct.zeta(Iy);
end
%% PLOTTING

[cols,nams]=flatcolors; colx = .5*cols{1}+.5*cols{21}; coly = cols{14}*.75+cols{24}*.25; col0 = cols{8}*.5+cols{4}*.5;
set_figsize([],12,12);
%Dipole moments
subaxis(2,1,1,'SV',0.004,'MT',0.025,'MR',.025,'MB',.08)
plot(minmax(elist),maxX.oscstrength(Nab)*[1,1],':','color',col0); hold on; 
px=plot(elist,maxX.oscstrength,'.-','color',colx,'LineWidth',.85,'MarkerSize',7); 
py=plot(elist,maxY.oscstrength,'.-','color',coly,'LineWidth',.85,'MarkerSize',7); hold off
set(gca,'Fontsize',8,'XTickLabel',{''},'LineWidth',.1)
ylim([0.34,1])
ylabel('Oscillator strength\ \ $|\langle \mathbf{r}\cdot \hat{\mathbf{n}}|\rho\rangle|^2$','Fontsize',10)
legend([px,py],{'$x$-polarized','$y$-polarized'},'Fontsize',10,'Orientation','Horizontal','Location','South'); legend boxoff

%Zeta
subaxis(2,1,2)
plot(minmax(elist),maxX.zeta(Nab)*[1,1],':','color',cols{8}*.5+cols{4}*.5); hold on; 
plot(elist,maxX.zeta,'.-','color',colx,'LineWidth',.85,'Markersize',7); 
plot(elist,maxY.zeta,'.-','color',coly,'LineWidth',.85,'Markersize',7); hold off
set(gca,'Fontsize',8,'LineWidth',.1)
ylabel('Eigenvalue\ \ $\zeta$','Fontsize',10)
xlabel('$\text{sgn}(a-b) \times \text{eccentricity}$','Fontsize',10)
ylim([0,7.15])

%Insets
inset = [1+5,Nab,numel(alist)-5]; xpos = [.09,.48,.84];
for ss = 1:3
    handaxes2 = axes('Position', [xpos(ss) .485 [0.2 0.135]*.6]);
    htp=patch(ellipse{inset(ss)}.bndry(:,1),ellipse{inset(ss)}.bndry(:,2),cols{19}*.9,'EdgeColor','none'); axis equal;
    xlim([-1,1]*1.07); ylim([-1,1]*1.07); axis off
end

break
%Save figure
path = 'C:\Data\Dropbox (MIT)\Documents\Notes\Maximum absorption in graphene nanostructures\img\matlabfrag\raw\';
name = 'DipoleResonanceEigenvalueParameters_Ellipse';
matlabfrag([path name])
fix_lines([path name '.eps'])

