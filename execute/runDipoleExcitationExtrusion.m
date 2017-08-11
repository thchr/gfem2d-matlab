function runDipoleExcitationExtrusion(bndshape,addstr,models,repincl,R,adivd,cut,Ncirc,plotresults)

addpath(genpath('../'));
ConstantsUnits0;

%%  ALLOW MULTITHREADING ON CLUSTER
fprintf('\n|----- MULTITHREADING SETUP -----|\n'); nopool = 1; %Don't open a multiworker pool
SetMatlabMultithreading;

%% DEFAULT INPUT
if ~exist('bndshape','var') || isempty(bndshape);
    bndshape = 'rectanglewedge'; %See posibilities under the switch statement for bndshape
end
if ~exist('addstr','var')
    addstr = [];
else
    addstr = ['_' addstr];
end

if ~exist('repincl','var') || isempty(repincl);
    repincl = [19,15]; %[15,13];
end
if ~exist('R','var') || isempty(R);
    R{1} = [1,0]; R{2} = [cosd(60),sind(60)];
end
if ~exist('adivd','var') || isempty(adivd);
    adivd = 2;
end
if ~exist('cut','var') || isempty(cut);
    cut = .401;
end
if ~exist('Ncirc','var') || isempty(Ncirc);
    Ncirc = 36;
end
if ~exist('plotresults','var') || isempty(plotresults);
    plotresults = 0;
end

%% ----- SETUP OF BOUNDARY AND MESHING -----

rep = repincl-cut*2;
repincl = bsxfun(@minus,repincl*2,mod(repincl,2));

%boundary of the ribbon unit cell, which consist of several lattice unit cells (first x-coords, then y-coords)
repcut = ceil(rep(1)*.65)+(1-cut)/adivd/2 + 0;
switch bndshape
    case 'rectangle'
        bnd = [zeros(1,2);
            [0,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)];
            rep(1)*R{1} + rep(2)*R{2};
            [rep(1)*R{1}(:,1) + rep(2)*R{2}(:,1),0]];
    case 'rectanglewedge'
        w=2-cut/2; h = 2*R{2}(2); repcut = ceil(rep(1)*.65)+(1-cut)/adivd/2 + 0;
        bnd = [zeros(1,2);
            [0,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)];
            [repcut*R{1}(:,1) - w/2,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)];
            [repcut*R{1}(:,1), rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2) - h];
            [repcut*R{1}(:,1) + w/2,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)];
            rep(1)*R{1} + rep(2)*R{2};
            [rep(1)*R{1}(:,1) + rep(2)*R{2}(:,1),0]];
end
H = 20; W=3+.5-cut/2; dw = 1-cut/2; dh = 1*R{2}(2);
bnd = [bnd; [repcut*R{1}(:,1) + W/2 + dw/2, 0]
            [repcut*R{1}(:,1) + W/2,-dh]; %add the "extrusion"
            [repcut*R{1}(:,1) + W/2, - H];
            [repcut*R{1}(:,1) - W/2,  - H];
            [repcut*R{1}(:,1) - W/2,-dh]; 
            [repcut*R{1}(:,1) - W/2 - dw/2,0]];
bnd = bsxfun(@minus,bnd, (rep(1)*R{1}+rep(2)*R{2})/2); %Re-center
bnd = unique(bnd,'rows','stable'); %Remove superfluous boundary points

d = bnd - circshift(bnd,[-1,0]); %Compute the boundary length of the structure under consideration
disp(['Boundary length = ' num2str(sum(sqrt(d(:,1).^2 + d(:,2).^2)),'%.1f')])

%% ----- DIPOLE PROPERTIES -----
r0 = [-4.5*R{1} + [0,rep(2)*R{2}(:,2)/2] + [-.5,.15],.125]; %position
polar = { [0,1,0] }; %polarization (y-polarized)

%% ----- AUXILIARY SQUARE MESH FOR PLOTTING -----
bndspc = .65; %boundary spacing

xstartend = [r0(1)-2.5,r0(1) + 12];
ystartend = [max(bnd(:,2)) - 3.5498,max(bnd(:,2))+1.25*bndspc];
xyfrac = ((max(bnd(:,2))+1.25*bndspc) - (R{2}(:,2)*2) ) / (xstartend(2) - xstartend(1));
x = linspace(xstartend(1),xstartend(2),400);
y = linspace(ystartend(1),ystartend(2),ceil(300*xyfrac));

[X,Y] = meshgrid(x,y);
%% ----- INCLUSIONS -----
%inclusion basic shape
theta = linspace(0,2*pi,Ncirc); theta=theta(1:end-1).'+theta(2)/2;
circ = [cos(theta),sin(theta)]/(2*adivd);

%--- regions where we reduce the number of points around the circle ---
%left upper corners
Ncircmod(1) = 20;
modregion{1}.x = [min(bnd(:,1))+[.75,.75],r0(1) - [3,3]];
modregion{1}.y = max(bnd(:,2))+[-4,-.75,-.75,-4];

%right upper corners
Ncircmod(end+1) = 20;
modregion{end+1}.x = [r0(1) + [12.5,12.5],max(bnd(:,1))-[.75,.75],];
modregion{end}.y = max(bnd(:,2))+[-4,-.75,-.75,-4];

%entire bottom region
Ncircmod(end+1) = 20;
modregion{end+1}.x = [min(bnd(:,1)),min(bnd(:,1)),max(bnd(:,1)),max(bnd(:,1))]+[1,1,-1,-1]*.75;
modregion{end}.y = [bnd(1,2),bnd(2,2),bnd(2,2),bnd(1,2)]+[.75,-4,-4,.75];

%bottom extrusion
Ncircmod(end+1) = 20;
modregion{end+1}.x = [bnd(end-1,1),bnd(end-1,1),bnd(end-3,1),bnd(end-3,1)] + [1,1,-1,-1]*.5;
modregion{end}.y = [bnd(end-2,2),bnd(end-1,2),bnd(end-1,2),bnd(end-2,2)]+ 0*[0,1,1,0]*dh + [.75,.5,.5,.75];

%left, right and bottom edges of holes
Ncircmod(end+1:end+3) = 28;
modregion{end+1}.x = [bnd(1,1),bnd(1,1),max(bnd(:,1)),max(bnd(:,1))]; %bottom
modregion{end}.y =    bnd(1,2)+[0,.75,.75,0];
modregion{end+1}.x = bnd(1,1) + [0,0,.75,.75]; %left
modregion{end}.y = [modregion{3}.y(1),modregion{1}.y(2),modregion{1}.y(2),modregion{3}.y(1)];
modregion{end+1}.x = max(bnd(:,1)) + [-.75,-.75,0,0]; %right
modregion{end}.y = [modregion{3}.y(1),modregion{2}.y(2),modregion{2}.y(2),modregion{3}.y(1)];

%define a (local) hdata density for each modregion
hdatafac = [[1,1,1]*6,6,[1,1,1]*4]; %"multiplication factor" for these regions
HeaviS = @(vals,val0,cut) (1 + erf(-(vals-val0)./cut))/2;
for cc = 1:numel(Ncircmod)
    hdatamod{cc} = @(x,y) HeaviS(-x,-modregion{cc}.x(1)-.75,1).*HeaviS(x,modregion{cc}.x(end)-.75,1).*...
        HeaviS(-y,-modregion{cc}.y(1)-.75,1).*HeaviS(y,modregion{cc}.y(2)-.75,1)*hdatafac(cc);
    
end

%---- create associated (origin-centered) circles for inclusions ----
for cc = 1:numel(Ncircmod)
    thetamod = linspace(0,2*pi,Ncircmod(cc)); thetamod=thetamod(1:end-1).'+thetamod(2)/2;
    circmod{cc} = [cos(thetamod),sin(thetamod)]/(2*adivd);
end

%add inclusions inside the unit cell
ii = 1;
for aa = ( (-10:repincl(1)+10) - (repincl(1)+1)/2 )
    for bb = ( (-20:repincl(2)) - (repincl(2)+1)/2 )
        center = aa*R{1} + bb*R{2};
        for cc = 1:numel(circmod)
            checkpos = inpolygon(center(:,1),center(:,2),modregion{cc}.x,modregion{cc}.y);
            if checkpos == 1
                circtype = cc; break
            else
                circtype = 0;
            end
        end
        if circtype == 0
            incl{ii} = bsxfun(@plus,circ, center);
        else
            incl{ii} = bsxfun(@plus,circmod{circtype}, center);
        end
        ii = ii + 1;
    end
end

%compute the (global) hdata function for the structure (mesh element size)
minincldist = 0.044819654451715;
hdata = @(x,y) (hdatamod{1}(x,y)+hdatamod{2}(x,y)+hdatamod{3}(x,y) + ...
    hdatamod{4}(x,y)+hdatamod{5}(x,y)+hdatamod{6}(x,y) + ...
    hdatamod{7}(x,y) + ...  %modified regions
    1) ... %elsewhere
    *minincldist*1.35;

if plotresults == 1
    %hdata function illustration
    set_figsize([],25,25);
    [Xhdata,Yhdata] = meshgrid(linspace(min(bnd(:,1)),max(bnd(:,1)),100),linspace(min(bnd(:,2)),max(bnd(:,2)),100));
    pcolor(Xhdata,Yhdata,hdata(Xhdata,Yhdata)); colormap(flipud(morgenstemning)); shading interp; 
    hold on;
    
    %limits of 'modregions'
    cols = flatcolors();
    plot(bnd([1:end,1],1),bnd([1:end,1],2),'-k'); hold on
    for cc = 1:numel(Ncircmod)
        plot(modregion{cc}.x([1:end,1]),modregion{cc}.y([1:end,1]),':k','color',cols{cc});
    end
    plot(xstartend([1,1,2,2,1]),ystartend([1,2,2,1,1]),'-','color',cols{14},'linewidth',1.75)
    plot(r0(1),r0(2),'p','MarkerSize',12.5,'MarkerFaceColor',cols{10}*.65+cols{6}*.35,'MarkerEdgeColor',cols{4})
    hold off; axis equal; drawnow
end

%% ----- MESH THE ACTUAL STRUCTURE -----
[p,t,nodes,cnct] = geomFancyInclusion(bnd,incl,hdata,plotresults);

%% ----- EXTERNAL POTENTIAL DUE TO DIPOLE -----

phiext = zeros(size(p,1),numel(polar)); %preallocate
for pp = 1:numel(polar)
    phiext(:,pp) = calcDipoleSource([p,zeros(size(p,1),1)],r0,polar{pp});
end

if strcmpi(plotresults,'yes') || all(plotresults == 1)
    [cols,nams]=flatcolors(); hold on;
    plot(r0(:,1),r0(:,2),'p','MarkerSize',12.5,'MarkerFaceColor',cols{10}*.65+cols{6}*.35,'MarkerEdgeColor',cols{4}); hold on;
    drawnow;
end

%% ----- EXTERNAL POTENTIAL ON AUXILIARY SQUARE MESH -----

phiextsq = zeros(numel(X),numel(polar)); %preallocate
for pp = 1:numel(polar)
    phiextsq(:,pp) = calcDipoleSource([X(:),Y(:),zeros(numel(X),1)],r0,polar{pp});
end

if strcmpi(plotresults,'yes') || all(plotresults == 1)
    hold on;
    plot([min(x),max(x),max(x),min(x),min(x)],[min(y),min(y),max(y),max(y),min(y)],'-','Color',cols{6});
    drawnow;
    hold off;
end

%% ----- MODEL(S) SETUP -----
if ~exist('models','var') || isempty(models)
    models = eval('models_nearzeroloss(4)'); %Near zero-loss at B = 4 T
elseif ischar(models)
    modelscall = models; clear models;
    if any(modelscall == '(') %In this case, the call is to a function, and a variable must be assigned
        models = eval(modelscall);
    else %In this case, the call is assumed to be to a script; variables automatically assigned
        eval(modelscall)
    end
elseif isnumeric(models) %Then we'll choose to interpret it as a B value in Tesla
    Bval = models; clear models
    ef_eV = 0.2; L = 400e-9;
    models = struct('ef_eV',ef_eV,'L',L,'gam_eV',1e-3 ,'ene_eV',60.3e-3,'B',Bval);
    addstr = [addstr '_B' num2str(Bval)];
end

%Print model type
fprintf('%g setup-models are considered:\n',numel(models))
modelsprint = struct2table(models(:)); disp(modelsprint); clear modelsprint
%% ----- PRINT THE SAVE NAME MODIFIER IF IT EXISTS -------
if ~isempty(addstr)
    fprintf('Savename modifier is: %s\n\n',addstr);
end
%% ----- CONSTRUCTING COULOMB MATRIX -----
fprintf('Estimate of necessary memory for matrix Coulomb matrix V: %g GB\n\n',size(p,1)^2*64/8*1e-9)

V = calcCoulomb(p,t);
fprintf('   Storing Coulomb matrix V requires %g GB memory\n\n',getfield(whos('V'),'bytes')/1e9)
%% ----- CONSTRUCTING REMAINING RELEVANT MATRICES AND SOLVING -----
%preallocate
results(size(models,1),size(models,2)) = struct('sigma',[],'sigma0',[],'rho',[],'phiind',[],'phiindsq',[]);

%setup parameters and conductivity functions
for mm=1:numel(models)
    fprintf('Model %g / %g [ omega = %.3f THz | L = %g nm | gamma = %g meV | B = %g T ] \n',...
        mm, numel(models), models(mm).ene_eV*eV2THz, models(mm).L*1e9, models(mm).gam_eV*1e3, models(mm).B );
    
    if ~isfield(models(mm),'f')
        results(mm).sigma0 = conducIntra(models(mm).ene_eV,models(mm).ef_eV,models(mm).gam_eV);
        results(mm).sigma = conducMagnetoClassical(models(mm).ene_eV,models(mm).ef_eV,models(mm).B,models(mm).gam_eV);
        f = results(mm).sigma/results(mm).sigma0;
    else
        results(mm).sigma0 = models(mm).sigma0;
        f = models(mm).f(p);  %Evaluate predefined function for 'f' on vertex points
    end
    models(mm).fvals = f;
    
    %weak form matrices
    [DR,DL] = calcDifferential(p,t,f);
    
    %setting up and solving system matrices with external potential
    omega = models(mm).ene_eV/hbar_eV;
    
    %Solve the matrix system
    ticslv = tic;
    results(mm).rho = ( (DL+1i*results(mm).sigma0/(4*pi*eps0*models(mm).L*omega)*DR*V) ) \ ...
                      ( -1i*results(mm).sigma0/(omega*models(mm).L^2)*DR*phiext );
    fprintf('   Matrix system solved in %.2f min\n\n',toc(ticslv)/60)
    
    results(mm).phiind = models(mm).L/(4*pi*eps0)*V*results(mm).rho;
end

%% REFINE THE RESULTS AROUND THE ZOOM-IN MESH ON A SQUARE MESH
clear V
Vsq = calcCoulombGeneralPoints(p,t,[X(:),Y(:)]); whosVsq = whos('Vsq');
fprintf('   Storing square-gridded Coulomb matrix Vsq requires %g GB memory\n\n',whosVsq.bytes/1e9)

for mm = 1:numel(models)
    results(mm).phiindsq = Vsq*results(mm).rho*models(mm).L/(4*pi*eps0);
end
clear Vsq
%% ----- SAVE DATA -----
savevars = {'results','models','phiext','phiextsq','p','t',...
    'nodes','cnct','x','y','r0','polar','repincl','Ncirc','bndshape'};
datapath = '../output/dipole/';
name = ['dipoleExciteLarge_' bndshape '_extrusion' addstr];
save([datapath name],savevars{:})
%% ----- PLOT RESULTS -----
if strcmpi(plotresults,'yes') || all(plotresults == 1)
    for mm = 1:numel(models)
        figure(mm)
        for pp = 1:numel(polar)
            subaxis(numel(polar),1,pp)
            v = results(mm).phiindsq(:,pp);
            
            set_figsize(4,25,25)
            %induced potential
            contourf(x*models(mm).L*1e6,y*models(mm).L*1e6,real(reshape(v,size(X))),101,'EdgeColor','none');
            set(gca,'YDir','normal')
            caxis(minmax(real(v))*.75); colormap(bluewhitered_mod)
            %outline of structure
            hold on
            for cc = 1:size(cnct,1)
                plot(nodes(cnct(cc,:),1)*models(mm).L*1e6,nodes(cnct(cc,:),2)*models(mm).L*1e6,'-','color',cols{4}*.5+cols{8}*.5,'LineWidth',.2)
            end
            %position of dipole source
            plot(r0(:,1)*models(mm).L*1e6,r0(:,2)*models(mm).L*1e6,'p','MarkerSize',12.5,'MarkerFaceColor',cols{10}*.65+cols{6}*.35,'MarkerEdgeColor',cols{4}*.5+cols{8}*.5);
            hold off
            axis on equal;
            xlim(minmax(x)*models(mm).L*1e6); ylim(minmax(y)*models(mm).L*1e6)
            set(gca,'Fontsize',8,'LineWidth',.2)
            xlabel('\itx\rm  [\mum]','Fontsize',10)
            ylabel('\ity\rm  [\mum]','Fontsize',10)
            drawnow;
        end
    end
end