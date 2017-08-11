clear all; clc; close all;
addpath(genpath('..'));

ConstantsUnits0;

cut = .401;
repincl =  [13,9];% [11,7];
Ncirc = 36;  %11
adivd = 2;
bndshape = 'rectanglewedge';
auxmesh = 'zoom';
usedefect = 'no';
plotresults = 'no';
collisiontest = 0;
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
    case 'rectangledefect'
        w=1-cut/2; h = 2*R{2}(2); repcut = ceil(rep(1)*.75)+(1-cut)/adivd/2 + 0;
        bnd = [zeros(1,2);
            [0,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)];
            [repcut*R{1}(:,1) - w/2,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)]; [repcut*R{1}(:,1) - w/2,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2) - h/50];
            [repcut*R{1}(:,1) - w/2,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2) - h];
            [repcut*R{1}(:,1) + w/2,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2) - h]; [repcut*R{1}(:,1) + w/2,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2) - h/50]
            [repcut*R{1}(:,1) + w/2,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)];
            rep(1)*R{1} + rep(2)*R{2};
            [rep(1)*R{1}(:,1) + rep(2)*R{2}(:,1),0]];
    case 'rectanglewedge'
        w=2-cut/2; h = 2*R{2}(2); repcut = ceil(rep(1)*.85)+(1-cut)/adivd/2 + 0;
        bnd = [zeros(1,2);
            [0,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)];
            [repcut*R{1}(:,1) - w/2,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)];
            [repcut*R{1}(:,1), rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2) - h];
            [repcut*R{1}(:,1) + w/2,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)];
            rep(1)*R{1} + rep(2)*R{2};
            [rep(1)*R{1}(:,1) + rep(2)*R{2}(:,1),0]];
    case 'rectangleasymwedge'
        w=4-cut/2; h = 4*R{2}(2); repcut = ceil(rep(1)*.85)+(1-cut)/adivd/2 + 0;
        bnd = [zeros(1,2);
            [0,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)];
            [repcut*R{1}(:,1) - w/2,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)];
            [repcut*R{1}(:,1) , rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2) - h];
            [repcut*R{1}(:,1) ,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)];
            rep(1)*R{1} + rep(2)*R{2};
            [rep(1)*R{1}(:,1) + rep(2)*R{2}(:,1),0]];
end
bnd = bsxfun(@minus,bnd, (rep(1)*R{1}+rep(2)*R{2})/2); %Re-center
bnd = unique(bnd,'rows','stable'); %Remove superfluous boundary points

%inclusion basic shape
theta = linspace(0,2*pi,Ncirc); theta=theta(1:end-1).'+theta(2)/2;
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
    w = 1; h = 2; Nd = 3;
    defect(:,1) = [linspace(w/2,-w/2,2).'; linspace(-w/2,-w/2,Nd).';  linspace(-w/2,w/2,2).';   linspace(w/2,w/2,Nd).' ];
    defect(:,2) = [linspace(h,h,2).'; linspace(h,-h,Nd).'; linspace(-h,-h,2).'; linspace(-h,h,Nd).'];
    defect = unique(defect,'rows','stable');
    defect = bsxfun(@plus,0.5*R{1} + [0,rep(2)*R{2}(:,2)/2],defect); %Move to boundary
    
    origdefect = defect;
    
    %test if any inclusions overlap with defect; if so, we have to handle that
    count = 1;
    for ii = 1:numel(incl)
        indefect = inpolygon(incl{ii}(:,1),incl{ii}(:,2),defect(:,1),defect(:,2));
        if any(indefect);
            rmvincl(count).cellind = ii;
            rmvincl(count).arrayind = find(~indefect);
            count = count + 1;
        end
    end
    for rr = 1:numel(rmvincl); %Remove points from defect list which lie within an inclusion
        definincl = inpolygon(defect(:,1),defect(:,2),incl{rmvincl(rr).cellind}(:,1),incl{rmvincl(rr).cellind}(:,2));
        defect = defect(~definincl,:);
    end
    if collisiontest == 1 %This entire thing is buggy under many circumstances, especially when it intersects in "non-convex hull"-like ways
        for rr = 1:numel(rmvincl) %Add the inclusions to the defect list
            defect = [defect; incl{rmvincl(rr).cellind}(rmvincl(rr).arrayind,:)];
        end
        
        unsort = defect(2:end,:); defect = defect(1,:);
        for dd = 1:size(unsort,1) %Sort the new defect list
            temp = unsort; conflict = 1;
            while conflict == 1 && ~all(isnan(temp(:))) %We don't want lines that go across the original defect
                [~,I] = min(pdist2(defect(dd,:),temp));
                midpoint = (defect(dd,:)+temp(I,:))/2;
                [IN,ON] = inpolygon(midpoint(:,1),midpoint(:,2),origdefect(:,1),origdefect(:,2));
                conflict = IN*~ON;
                if conflict == 1
                    for rr = 1:numel(rmvincl)
                        midpoint = (defect(dd,:)+temp(I,:))/2;
                        [IN,ON] = inpolygon(midpoint(:,1),midpoint(:,2),incl{rmvincl(rr).cellind}(:,1),incl{rmvincl(rr).cellind}(:,2));
                        conflict = IN*~ON;
                        if conflict == 0; break; end
                    end
                end
                hold off; axis equal; drawnow; pause
                temp(I,:) = NaN;
            end
            defect(dd+1,:) = unsort(I,:);
            unsort(I,:) = [];
        end
    end
    %Remove the inclusions that overlapped with the original defect
    for rr = numel(rmvincl):-1:1;  incl(rmvincl(rr).cellind) = [];  end
    
    %Finally, add the new composite defect to the inclusion list
    incl{end+1} = defect;
    
    %Text for saving
    defecttxt = 'WithDefect';
else
    defecttxt = 'NoDefect';
end

%mesh the structure
[p,t,nodes,cnct] = geomFancyInclusion(bnd,incl,@(x,y) 20*ones(size(x)),0);

%% ----- EXTERNAL POTENTIAL -----
r0 = [-1.5*R{1} + [0,rep(2)*R{2}(:,2)/2] + [0,.15],.125]; %position
polar = {[1,0,0],[0,1,0],[0,0,1]};
phiext = zeros(size(p,1),numel(polar)); %preallocate
for pp = 1:numel(polar)
    phiext(:,pp) = calcDipoleSource([p,zeros(size(p,1),1)],r0,polar{pp});
end

if strcmpi(plotresults,'yes') || all(plotresults == 1)
    [cols,nams]=flatcolors();
    plot(r0(:,1),r0(:,2),'p','MarkerSize',12.5,'MarkerFaceColor',cols{10}*.65+cols{6}*.35,'MarkerEdgeColor',cols{4}); hold on;
    trisurf(t,p(:,1),p(:,2),phiext,'FaceAlpha',.5,'EdgeColor',cols{4},'LineWidth',.1);
    view(2); colormap(bluewhitered_mod); caxis(minmax(phiext)/20); axis equal;  hold off;
    drawnow;
end

%% ----- AUXILIARY SQUARE MESH -----
bndspc = .65; %boundary spacing
switch auxmesh
    case 'full'
        x = linspace(min(p(:,1))-bndspc,max(p(:,1)+bndspc),200); y = linspace(min(p(:,2))-bndspc,max(p(:,2)+bndspc),250);
    case 'zoom'
        x = linspace(r0(1)-4*bndspc,max(p(:,1)+bndspc),200); y = linspace(0,max(p(:,2)+bndspc),200);
end
[X,Y] = meshgrid(x,y);
phiextsq = zeros(numel(X),numel(polar)); %preallocate
for pp = 1:numel(polar)
    phiextsq(:,pp) = calcDipoleSource([X(:),Y(:),zeros(numel(X),1)],r0,polar{pp});
end
%% ----- CONSTRUCTING COULOMB MATRIX -----
V = calcCoulomb(p,t);
%Vsq = calcCoulombGeneralPoints(p,t,[X(:),Y(:)]);

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
%LHS =  (DL+1i*sigma0/(4*pi*eps0*L*omega)*DR*V);
%RHS = -1i*sigma0/(omega*L^2)*DR*phiext;
ticslv = tic; 
rho = ( (DL+1i*sigma0/(4*pi*eps0*L*omega)*DR*V) ) \ ( -1i*sigma0/(omega*L^2)*DR*phiext );
fprintf('\nMatrix system solve in %.2f min\n\n',toc(ticslv)/60)
%rho = LHS\RHS;
phiind = L/(4*pi*eps0)*V*rho;


clear V
Vsq = calcCoulombGeneralPoints(p,t,[X(:),Y(:)]);
phiindsq = Vsq*rho*L/(4*pi*eps0);
%% ----- SAVE DATA -----
savevars = {'phiind','rho','phiext','phiindsq','phiextsq','p','t',...
    'nodes','cnct','ene_eV','ef_eV','gam_eV','B','L','x','y','r0','polar'};
datapath = '../output/dipole/';
name = ['dipExciteOneWay_gam' strrep(num2str(gam_eV*1e3),'.','p') 'meV_B' num2str(B) 'T_' bndshape '_' auxmesh];
save([datapath name],savevars{:})
%% ----- PLOT RESULTS -----
if strcmpi(plotresults,'yes') || all(plotresults == 1)
    for pp = 1:numel(polar)
        subaxis(numel(polar),1,pp)
        v = phiindsq(:,pp);
        
        set_figsize(4,25,25)
        %induced potential
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
    end
    export_fig(['../figures/dipole/' savename],'-pdf')
end