clear all; close all; clc;

ConstantsUnits0;
cols = flatcolors; darkcol = cols{4};

%% MESH
geom = 'disk';
if strcmpi(geom,'disk')
    mesh = geomEllipse([1,1],[],.1,1);
elseif strcmpi(geom,'ellipse')
    mesh = geomEllipse([1,.75],[],.1,1);
elseif strcmpi(geom,'trilattice') %THIS SCENARIO DOES NOT MAKE SENSE TO TEST THE EDGE STATES BECAUSE THEY REQUIRE NONZERO INCIDENT WAVE VECTOR (AND PLANE WAVE IS NECESSARILY NORMAL INCIDENCE = ZERO WAVE VECTOR)
    %defining parameters
    R{1} = [1,0]; R{2} = [cosd(60),sind(60)];
    repincl = [3,7];
    cut = .401;
    %boundary and repetition index
    rep = repincl-cut*2;
    repincl = bsxfun(@minus,repincl.*4,mod(repincl,2));
    bnd = [zeros(1,2); [0,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)]; rep(1)*R{1} + rep(2)*R{2}; [rep(1)*R{1}(:,1) + rep(2)*R{2}(:,1),0]];
    bnd = bsxfun(@minus,bnd, (rep(1)*R{1}+rep(2)*R{2})/2); %Re-center
    %inclusion basic shape
    adivd = 2;   Ncirc = 16;
    theta = linspace(0,2*pi,Ncirc); theta=theta(1:end-1).'+theta(2)/2;
    circ = [cos(theta),sin(theta)]/(2*adivd);
    %add inclusions inside the unit cell
    ii = 1; clear incl
    for aa = ( (1:repincl(1)) - (repincl(1)+1)/2 )
        for bb = ( (1:repincl(2)) - (repincl(2)+1)/2 )
            incl{ii} = bsxfun(@plus,circ, aa*R{1} + bb*R{2});
            ii = ii + 1;
        end
    end
    %mesh the structure
    [mesh.p,mesh.t,nodes,cnct] = geomFancyInclusion(bnd,incl,@(x,y) ones(size(x))*.2,1);
end
[alpha,eigstruct] = calcPolarizability(mesh,20);
%% GRAPHENE
ef = .2;
kbt =kb_eV*300;
gam = 1e-3;
ene = linspace(0.001,.1,250); ene0 = linspace(0.001,.1,2500);
sigma = conducIntra(ene,ef,gam); sigma0 = conducIntra(ene0,ef,gam);

B = 4;
[~,sigmaxx,sigmaxy] = conducMagnetoClassical(ene,ef,B,gam);
f(1,1,:) = sigmaxx./sigma;
f(1,2,:) = sigmaxy./sigma;
f(2,1,:) = -sigmaxy./sigma;
f(2,2,:) = sigmaxx./sigma;
omegac = e*B*vf^2/(ef*ev2jo)*hbar_eV;
%% RESPONSE
L = 200e-9;
polar = [1,0;0,1;[1,1i]/sqrt(2);[1,-1i]/sqrt(2)]; %[x-linear, y-linear, left-circular, right-circular]
zeta0 = 2i*eps0*ene0/hbar_eV*L./sigma0;

Q_abs0 = (repmat(ene0(:),1,2)/hbar_eV/c).*imag([alpha.x(zeta0,L);alpha.y(zeta0,L)].')/(sum(meshArea(mesh.p,mesh.t))*L^2);
Q_abs = calcAbsorption(mesh.p,mesh.t,ene/hbar_eV,L,polar,sigma,f);

%% PLOT
set_figsize(2,20,15)
plot([1,1]*omegac,[0,2],':','color',cols{8}); hold on
plot(ene0,Q_abs0(:,1),'color',cols{4},'LineWidth',2)
plot(ene0,Q_abs0(:,2),'color',cols{4}*.5+cols{8}*.5,'LineWidth',2)
hold on
plot(ene,Q_abs(:,1),'-.','color',cols{1},'LineWidth',2) %x
plot(ene,Q_abs(:,2),':','color',cols{2},'LineWidth',2)  %y
plot(ene,Q_abs(:,3),':','color',cols{3},'LineWidth',2)  %left
plot(ene,Q_abs(:,4),':','color',cols{5},'LineWidth',2)  %right
hold off