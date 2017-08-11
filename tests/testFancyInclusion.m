clear all; clc; close all;
ConstantsUnits0; [cols,nams]=flatcolors();

cut = .4;
repincl = [11,9];
adivd = 2;
bndshape = 'rectangle';
auxmesh = 'zoom'; 
usedefect = 'yes';
R{1} = [1,0]; R{2} = [cosd(60),sind(60)];

rep = repincl-cut*2;
repincl = bsxfun(@minus,repincl*2,mod(repincl,2));

%% ----- SETUP OF BOUNDARY AND MESHING -----
%boundary of the ribbon unit cell, which consist of several lattice unit cells (first x-coords, then y-coords)
switch bndshape
    case 'trapezoid'
        bnd = [zeros(1,2); rep(2)*R{2}; rep(1)*R{1}+rep(2)*R{2}; rep(1)*R{1}];
    case 'rectangle'
        bnd = [zeros(1,2); [0,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)]; rep(1)*R{1} + rep(2)*R{2}; [rep(1)*R{1}(:,1) + rep(2)*R{2}(:,1),0]];
end
bnd = bsxfun(@minus,bnd, (rep(1)*R{1}+rep(2)*R{2})/2); %Re-center
bnd = unique(bnd,'rows','stable'); %Remove superfluous boundary points

%inclusion basic shape
theta = linspace(0,2*pi,14); theta=theta(1:end-1).'+theta(2)/2;
circ = [cos(theta),sin(theta)]/(2*adivd);

%add inclusions inside the unit cell
ii = 1;
for aa = ( (1:repincl(1)) - (repincl(1)+1)/2 )
    for bb = ( (1:repincl(2)) - (repincl(2)+1)/2 )
        incl{ii} = bsxfun(@plus,circ, aa*R{1} + bb*R{2});
        ii = ii + 1;
    end
end

%make a defect
if strcmpi(usedefect,'yes') || all(usedefect == 1)
    w = .15; h = .5; Nd = 3;
    defect(:,1) = [linspace(w/2,-w/2,2).'; linspace(-w/2,-w/2,Nd).';  linspace(-w/2,w/2,2).';   linspace(w/2,w/2,Nd).' ];
    defect(:,2) = [linspace(h,h,2).'; linspace(h,-h,Nd).'; linspace(-h,-h,2).'; linspace(-h,h,Nd).'];
    defect = bsxfun(@plus,0.5*R{1} + [0,rep(2)*R{2}(:,2)/2],defect); %Move to boundary
    incl{end+1} = defect;
end

%mesh the structure
[p,t,nodes,cnct] = geomFancyInclusion(bnd,incl,@(x,y) 15*ones(size(x)),'mesh');

%% ----- EXTERNAL POTENTIAL -----
r0 = [-2.5*R{1} + [0,rep(2)*R{2}(:,2)/2] - [0,.15],.125]; %position
polar = [1,0,0];
phiext = calcDipoleSource([p,zeros(size(p,1),1)],r0,polar);

hold on; plot(r0(:,1),r0(:,2),'p','MarkerSize',12.5,'MarkerFaceColor',cols{10}*.65+cols{6}*.35,'MarkerEdgeColor',cols{4}); hold off; drawnow

%% ----- AUXILIARY SQUARE MESH -----
bndspc = .65; %boundary spacing
switch auxmesh 
    case 'full'
        x = linspace(min(p(:,1))-bndspc,max(p(:,1)+bndspc),200); y = linspace(min(p(:,2))-bndspc,max(p(:,2)+bndspc),125);
    case 'zoom'
        x = linspace(r0(1)-4*bndspc,max(p(:,1)+bndspc),200); y = linspace(0,max(p(:,2)+bndspc),100);
end
[X,Y] = meshgrid(x,y);
phiextsq = calcDipoleSource([X(:),Y(:),zeros(numel(X),1)],r0,polar);

%% ----- CONSTRUCTING COULOMB MATRIX -----
V = calcCoulomb(p,t);
Vsq = calcCoulombGeneralPoints(p,t,[X(:),Y(:)]);

%% ----- CONSTRUCTING REMAINING RELEVANT MATRICES AND SOLVING -----
%setup parameters and conductivity functions
ef_eV = 0.2; gam_eV = 1e-3; ene_eV = 60e-3; B = 8; L = 400e-9;
sigma0 = conducLRA(ene_eV,ef_eV,gam_eV);
sigma = conducMagnetoClassical(ene_eV,ef_eV,B,gam_eV);
f = sigma/sigma0; %f=eye(2);

%weak form matrices
[DR,DL] = calcDifferential(p,t,f);

%setting up and solving system matrices with external potential
omega = ene_eV/hbar_eV;
LHS =  (DL+1i*sigma0/(4*pi*eps0*L*omega)*DR*V);
RHS = -1i*sigma0/(omega*L^2)*DR*phiext;
rho = LHS\RHS; phiind = V*rho*L/(4*pi*eps0);
phiindsq = Vsq*rho*L/(4*pi*eps0);

%% ----- PLOT RESULTS -----
v = phiindsq; 
set_figsize(4,25,25)

%induced potential
%imagesc(x,y,real(reshape(v,size(X))));
contourf(x*L*1e6,y*L*1e6,real(reshape(v,size(X))),101,'EdgeColor','none');
set(gca,'YDir','normal')
caxis(minmax(real(v))*.75); colormap(bluewhitered_mod)
%outline of structure
hold on
for cc = 1:size(cnct,1)
    plot(nodes(cnct(cc,:),1)*L*1e6,nodes(cnct(cc,:),2)*L*1e6,'-','color',cols{4}*.5+cols{8}*.5,'LineWidth',.2)
end
drawnow;
%position of dipole source
plot(r0(:,1)*L*1e6,r0(:,2)*L*1e6,'p','MarkerSize',12.5,'MarkerFaceColor',cols{10}*.65+cols{6}*.35,'MarkerEdgeColor',cols{4}*.5+cols{8}*.5); 
hold off
axis on equal; 
xlim(minmax(x)*L*1e6); ylim(minmax(y)*L*1e6)
set(gca,'Fontsize',8,'LineWidth',.2)
xlabel('\itx\rm  [\mum]','Fontsize',10)
ylabel('\ity\rm  [\mum]','Fontsize',10)

%export_fig('dipoleExcitedOneWay_WithDefect_VeryLowLoss','-pdf')