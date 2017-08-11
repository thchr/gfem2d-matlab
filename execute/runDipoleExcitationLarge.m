function runDipoleExcitationLarge(bndshape,auxmesh,addstr,models,repincl,R,adivd,cut,Ncirc,plotresults)

addpath(genpath('../'));
ConstantsUnits0;

%%  ALLOW MULTITHREADING ON CLUSTER
fprintf('\n|----- MULTITHREADING SETUP -----|\n'); nopool = 1; %Don't open a multiworker pool
SetMatlabMultithreading;

%% DEFAULT INPUT
if ~exist('bndshape','var') || isempty(bndshape);
    bndshape = 'rectangle'; %See posibilities under the switch statement for bndshape
end
if ~exist('auxmesh','var') || isempty(auxmesh);
    auxmesh = 'oneedgezoom';
end
if ~exist('addstr','var')
    addstr = [];
else
    addstr = ['_' addstr];
end

if ~exist('repincl','var') || isempty(repincl);
    repincl = [27,21]; %[15,13];
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
        w=2-cut/2; h = 2*R{2}(2); repcut = ceil(rep(1)*.65)+(1-cut)/adivd/2 + 0;
        bnd = [zeros(1,2);
            [0,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)];
            [repcut*R{1}(:,1) - w/2,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)];
            [repcut*R{1}(:,1), rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2) - h];
            [repcut*R{1}(:,1) + w/2,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)];
            rep(1)*R{1} + rep(2)*R{2};
            [rep(1)*R{1}(:,1) + rep(2)*R{2}(:,1),0]];
    case 'rectangleasymwedge'
        w=4-cut/2; h = 4*R{2}(2); repcut = ceil(rep(1)*.65)+(1-cut)/adivd/2 + 0;
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

d = bnd - circshift(bnd,[-1,0]); %Compute the boundary length of the structure under consideration
disp(['Boundary length = ' num2str(sum(sqrt(d(:,1).^2 + d(:,2).^2)),'%.1f')])

%inclusion basic shape
theta = linspace(0,2*pi,Ncirc); theta=theta(1:end-1).'+theta(2)/2;
circ = [cos(theta),sin(theta)]/(2*adivd);
%--- regions where we reduce the number of points around the circle ---
%right and left upper corners
Ncircrough1 = 16; 
thetarough1 = linspace(0,2*pi,Ncircrough1); thetarough1=thetarough1(1:end-1).'+thetarough1(2)/2;
circrough1 = [cos(thetarough1),sin(thetarough1)]/(2*adivd);
roughregion1.x = {max(bnd(:,1))-[.5,.5,11,11], min(bnd(:,1))+[11,11,.5,.5]}; 
roughregion1.y = {max(bnd(:,2))+[-4,-.75,-.75,-4], max(bnd(:,2))+[-4,-.75,-.75,-4]};
%entire bottom region
Ncircrough2 = 16; 
thetarough2 = linspace(0,2*pi,Ncircrough2); thetarough2=thetarough2(1:end-1).'+thetarough2(2)/2;
circrough2 = [cos(thetarough2),sin(thetarough2)]/(2*adivd);
roughregion2.x = [min(bnd(:,1)),max(bnd(:,1)),max(bnd(:,1)),min(bnd(:,1))]-[-1,1,1,-1]*.5; 
roughregion2.y = [min(bnd(:,2)),min(bnd(:,2)),max(bnd(:,2)),max(bnd(:,2))]-[-.75,-.75,4,4];
%left, right and bottom edges of holes
Ncircrough3 = 18; 
thetarough3 = linspace(0,2*pi,Ncircrough3); thetarough3=thetarough3(1:end-1).'+thetarough3(2)/2;
circrough3 = [cos(thetarough3),sin(thetarough3)]/(2*adivd);
roughregion3.x = {[min(bnd(:,1)),max(bnd(:,1)),max(bnd(:,1)),min(bnd(:,1))],... %bottom
                  min(bnd(:,1)) + [0,.75,.75,0],... %left
                  max(bnd(:,1)) + [-.75,0,0,-.75]}; %right
roughregion3.y = {min(bnd(:,2))+[0,0,.75,.75], roughregion2.y+[0,0,1.5,1.5], roughregion2.y+[0,0,1.5,1.5]}; %bottom, left, and right

%add inclusions inside the unit cell
ii = 1;
for aa = ( (1:repincl(1)) - (repincl(1)+1)/2 )
    for bb = ( (1:repincl(2)) - (repincl(2)+1)/2 )
        center = aa*R{1} + bb*R{2}; 
        if inpolygon(center(:,1),center(:,2),roughregion2.x,roughregion2.y)
            incl{ii} = bsxfun(@plus,circrough2, center);
        elseif inpolygon(center(:,1),center(:,2),roughregion1.x{1},roughregion1.y{1}) || ...
               inpolygon(center(:,1),center(:,2),roughregion1.x{2},roughregion1.y{2}) 
            incl{ii} = bsxfun(@plus,circrough1, center);
        elseif inpolygon(center(:,1),center(:,2),roughregion3.x{1},roughregion3.y{1}) || ...
               inpolygon(center(:,1),center(:,2),roughregion3.x{2},roughregion3.y{2}) || ...
               inpolygon(center(:,1),center(:,2),roughregion3.x{3},roughregion3.y{3}) 
            incl{ii} = bsxfun(@plus,circrough3, center);
        else
            incl{ii} = bsxfun(@plus,circ, center);
        end
        ii = ii + 1;
    end
end

%hdata function object (local mesh element size)
minincldist = 0.044819654451715;
HeaviS = @(vals,val0,cut) (1 + erf(-(vals-val0)./cut))/2; 
hdata = @(x,y) ((HeaviS(y,max(bnd(:,2))-6,2).*HeaviS(-y,-min(bnd(:,2))-.5,1).* ...
                 HeaviS(-x,-min(bnd(:,1))-1.25,1).*HeaviS(x,max(bnd(:,1))-1.25,1) + ... %~bottom half
                 HeaviS(y,max(bnd(:,2))-1.25,1).*HeaviS(-y,-(max(bnd(:,2))-5.5),1).*...
                 HeaviS(-x,-min(bnd(:,1))-1.25,1).*HeaviS(x,min(bnd(:,1))+9.5,1) + ... %left top
                 HeaviS(y,max(bnd(:,2))-1.25,1).*HeaviS(-y,-(max(bnd(:,2))-5.5),1).*...
                 HeaviS(-x,-(max(bnd(:,1))-8.75),1).*HeaviS(x,(max(bnd(:,1))-1.25),1)) ... %right top
                 *8 ... %"multiplication factor" for these regions
                +1) ... %elsewhere
                *minincldist*1.65;
% [X,Y] = meshgrid(linspace(min(bnd(:,1)),max(bnd(:,1)),100),linspace(min(bnd(:,2)),max(bnd(:,2)),100));
% surf(X,Y,hdata(X,Y)); colormap(hot)
% axis equal;
% hold on; plot(roughregion.x([1:end,1]),roughregion.y([1:end,1]),'-g'); drawnow
% return

%mesh the structure
[p,t,nodes,cnct] = geomFancyInclusion(bnd,incl,hdata,plotresults);

%% ----- EXTERNAL POTENTIAL -----
r0 = [-4.5*R{1} + [0,rep(2)*R{2}(:,2)/2] + [0,.15],.125]; %position
polar = { [0,1,0] }; %y-polarized

phiext = zeros(size(p,1),numel(polar)); %preallocate
for pp = 1:numel(polar)
    phiext(:,pp) = calcDipoleSource([p,zeros(size(p,1),1)],r0,polar{pp});
end

if strcmpi(plotresults,'yes') || all(plotresults == 1)
    [cols,nams]=flatcolors(); hold on; 
    plot(r0(:,1),r0(:,2),'p','MarkerSize',12.5,'MarkerFaceColor',cols{10}*.65+cols{6}*.35,'MarkerEdgeColor',cols{4}); hold on;
    drawnow;
end

%% ----- AUXILIARY SQUARE MESH -----
bndspc = .65; %boundary spacing
switch auxmesh
    case 'full'
        x = linspace(min(p(:,1))-bndspc,max(p(:,1)+bndspc),200); y = linspace(min(p(:,2))-bndspc,max(p(:,2)+bndspc),250);
    case 'zoom'
        x = linspace(r0(1)-4*bndspc,max(p(:,1)+bndspc),200); y = linspace(0,max(p(:,2)+bndspc),200);
    case 'oneedgezoom'
        xstartend = [r0(1)-3.85*bndspc,max(p(:,1))-4.075*bndspc]-[0,8];
        xyfrac = ((max(p(:,2))+1.25*bndspc) - (R{2}(:,2)*2) ) / (xstartend(2) - xstartend(1));
        x = linspace(xstartend(1),xstartend(2),400); 
        y = linspace(max(p(:,2)) - 3.5498,max(p(:,2))+1.25*bndspc,ceil(300*xyfrac));
end
[X,Y] = meshgrid(x,y);
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
if numel(models) > 1
    error('Because of memory constraints, we later overwrite the Coulomb matrix - accordingly, only one model is allowed\n')
end

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
    
    %We have to overwrite V because we run out of memory otherwise (this
    %sets up the LHS of the matrix system) - this means we CANNOT do
    %multiple models
    ticslhs = tic;
    V = ( (DL+1i*results(mm).sigma0/(4*pi*eps0*models(mm).L*omega)*DR*V) );
    fprintf('   LHS of matrix system computed in %.2f min\n',toc(ticslhs)/60)
    fprintf('      (requires %g GB memory)\n\n',getfield(whos('V'),'bytes')/1e9)
    
    %Solve the matrix system
    ticslv = tic; 
    results(mm).rho = V \ ( -1i*results(mm).sigma0/(omega*models(mm).L^2)*DR*phiext );
    fprintf('   Matrix system solved in %.2f min\n\n',toc(ticslv)/60)
    
    %results(mm).phiind = models(mm).L/(4*pi*eps0)*V*results(mm).rho;
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
name = ['dipoleExciteLarge_' bndshape '_' auxmesh addstr];
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