clc; clear all; close all;
%Beware, there seems to be a few errors hidden through (e.g. realted to
%prefactors and size-scaling). Noted on May 24, 2016.


addpath(genpath('..'));
ConstantsUnits0; [cols,nams] = flatcolors;
ncirc = 25;                                          %is a 1/4 radius disk
theta = linspace(0,2*pi,ncirc+1); theta = theta(1:end-1);
inclusion = .25*[cos(theta);sin(theta)].';
rep = [7,6];
[~,origmesh] = geomRibbonInclusion('triangular',rep,inclusion,20,@(x,y) 5*ones(size(x)),'origmesh');
p = origmesh.p; t = origmesh.t;
%%
V = calcCoulomb(p,t);
%%
x = linspace(min(p(:,1))-.5,max(p(:,1)+.5),250); y = linspace(min(p(:,2))-.5,max(p(:,2)+.5),250);
[X,Y] = meshgrid(x,y);
Vsq = calcCoulombGeneralPoints(p,t,[X(:),Y(:)]);
%%
set_figsize(2,20,20)
r0 = [.8,max(p(:,2))-.05,.125];
[DR,DL] = calcDifferential(p,t);
phiext = calcDipoleSource([p,zeros(size(p,1),1)],r0,[0,0,1]);
phiextsq = calcDipoleSource([X(:),Y(:),zeros(numel(X),1)],r0,[0,0,1]);
%phiext = -p(:,1);
ef_eV = .2; L = 100e-9; gam_eV = 6e-3; B = 10; eneB_eV = hbar_eV*e*B*vf^2/(ef_eV*ev2jo);
omega = linspace(.1,1.25,25)*ef_eV/hbar_eV;
sigma = @(om) conducLRA(om*hbar_eV,.2,gam_eV);
sigma1 = @(om) conducMagnetoClassical(om*hbar_eV,ef_eV,20,gam_eV);
f = @(om) sigma1(om)/sigma(om);
ticloop = tic;
for oo = 1:numel(omega)
    
    [~,DR] = evalc('calcDifferential(p,t,f(omega(oo)))'); %Run without output to command line
    LHS =  (DL+1i*sigma(omega(oo))/(4*pi*eps0*L*omega(oo))*DR*V);
    RHS = -1i*sigma(omega(oo))/(omega(oo)*L)*DR*phiext;
    rho = LHS\RHS; phi = V*rho/(4*pi*eps0*L); phiind(:,oo) = phi-phiext;
    phiindsq(:,oo) = Vsq*rho/(4*pi*eps0*L) - phiextsq;
    
    if mod(oo,1) == 0;
        fprintf('   %g/%g frequencies calculated (%.1f min/%.1f min)\n',oo,numel(omega),toc(ticloop)/60,toc(ticloop)/60*numel(omega)/oo);
    end
    
end

%%
R1 = [1,0]; R2 = [cosd(60),sind(60)];
bnd = [R1(:,1)*minmax((1:rep(1)) - (rep(1)+1)/2) + R2(:,1)*(rep(2)-1)/2 + ([-1,1]*R1(:,1)+[1,1]*R2(:,1))/2, fliplr( R1(:,1)*minmax((1:rep(1)) - (rep(1)+1)/2) - R2(:,1)*(rep(2)-1)/2 + ([-1,1]*R1(:,1)-[1,1]*R2(:,1))/2 );
    R1(:,2)*minmax((1:rep(1)) - (rep(1)+1)/2) + R2(:,2)*(rep(2)-1)/2 + ([-1,1]*R1(:,2)+[1,1]*R2(:,2))/2, fliplr( R1(:,2)*minmax((1:rep(1)) - (rep(1)+1)/2) - R2(:,2)*(rep(2)-1)/2 + ([-1,1]*R1(:,2)-[1,1]*R2(:,2))/2 )];
set_figsize(2,38,25)
for oo = 1:numel(omega);
    clims = .75*[-1,1]*max(abs(real(phiindsq(:,oo))));
    subaxis(5,5,oo,'M',.01,'MT',.0275,'SV',.0275,'SH',.0025)
    pcolor(X,Y,real(reshape(phiindsq(:,oo),size(X))));  hold on; shading interp;
    caxis(clims);
    plot3(r0(1),r0(2),clims(2)*1.1,'*','Color',cols{3}); 
    colormap('bluewhitered_mod');
    %colormap('blueblackred')
    
    
    %Plot inclusions
    for rr = (1:rep(1)) - (rep(1)+1)/2;
        for uu=(1:rep(2)) - (rep(2)+1)/2;
            if rr==0; lnst='-'; else lnst='-'; end %Linestyle
            plot(inclusion([1:end,1],1) + rr*R1(:,1) + uu*R2(:,1),...
                inclusion([1:end,1],2) + rr*R1(:,2) + uu*R2(:,2),...
                lnst,'Color',cols{8},'linewidth',.2); hold on;
        end
    end
    %Plot ribbon boundaries
    plot(bnd(1,[1:end,1]),bnd(2,[1:end,1]),'-','color',cols{8},'linewidth',.2)
    hold off;
    
    view(2);
    axis equal off;
    text(min(p(:,1)),max(p(:,2))*1.025,['\omega/\omega_c = ' num2str(omega(oo)*hbar_eV/eneB_eV,2)],...
        'Fontsize',9,'Color',cols{4},'VerticalAlignment','top')
    drawnow; pause(.01)
end

%%
zoom.x = [-2,5.25]; zoom.y = [1,2.85];
xZ = linspace(zoom.x(1),zoom.x(2), round(250*diff(zoom.x)/diff(zoom.y)) ); 
yZ = linspace(zoom.y(1),zoom.y(2), round(250*diff(zoom.y)/diff(zoom.x)) );
[XZ,YZ] = meshgrid(xZ,yZ);

VZ = calcCoulombGeneralPoints(p,t,[XZ(:),YZ(:)]);
phiextZ = calcDipoleSource([XZ(:),YZ(:),zeros(numel(XZ),1)],r0,[0,0,1]);
%%
omegaZ = eneB_eV*4.5/hbar_eV;
%--------------
[~,DR] = evalc('calcDifferential(p,t,f(omegaZ))'); %Run without output to command line
LHS =  (DL+1i*sigma(omegaZ)/(4*pi*eps0*L*omegaZ)*DR*V);
RHS = -1i*sigma(omegaZ)/(omegaZ*L)*DR*phiext;
rhoZ = LHS\RHS; phi = V*rhoZ/(4*pi*eps0*L);
phiindsqZ = VZ*rhoZ/(4*pi*eps0*L) - phiextZ;
%--------------
%%
set_figsize(3,25,25)
colg = cols{4}*.25 + cols{8}*.75;
clims = [-1,1]*max(abs(real(phiindsqZ)));
subaxis(1,1,1,'M',.1,'MT',.0275,'SV',.0275,'SH',.0025)
%pcolor(XZ,YZ,real(reshape(phiindsqZ,size(XZ))));  hold on; shading interp;
%imagesc(xZ,yZ,real(reshape(phiindsqZ,size(XZ))));  hold on; set(gca,'YDir','normal')
contourf(XZ,YZ,real(reshape(phiindsqZ,size(XZ))),linspace(clims(1),clims(2),150),'LineStyle','none'); hold on

caxis(.75*clims);

colormap('bluewhitered_mod');
%colormap('blueblackred')
axis equal off;

%Plot inclusions
for rr = (1:rep(1)) - (rep(1)+1)/2;
    for uu=(1:rep(2)) - (rep(2)+1)/2;
        if rr==0; lnst='-'; else lnst='-'; end %Linestyle
        plot(inclusion([1:end,1],1) + rr*R1(:,1) + uu*R2(:,1),...
            inclusion([1:end,1],2) + rr*R1(:,2) + uu*R2(:,2),...
            lnst,'Color',colg,'linewidth',.3); hold on;
    end
end
%Plot ribbon boundaries
plot(bnd(1,[1:end,1]),bnd(2,[1:end,1]),'-','color',colg,'linewidth',.3)
plot3(r0(1),r0(2),1,'p','Color',colg,'MarkerFaceColor',cols{3},'MarkerSize',11,'LineWidth',.1); 
hold off;

drawnow
xlim(zoom.x)
ylim(zoom.y)

%export_fig('DipoleExcitation_OneWayEdgePropagation','-pdf');


