function runDipoleTimeEvolution(bndshape,addstr,models,excitation,repincl,R,adivd,cut,Ncirc,plotresults)

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

if ~exist('excitation','var') || isempty(excitation)
    excitation.center_eV = 60.3e-3;
    excitation.width_eV = 5e-3; %~ 5 THz
    excitation.profile = 'pulse';
end

if ~exist('repincl','var') || isempty(repincl);
    repincl = [15,13];
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

bnd = bsxfun(@minus,bnd, (rep(1)*R{1}+rep(2)*R{2})/2); %Re-center
bnd = unique(bnd,'rows','stable'); %Remove superfluous boundary points

d = bnd - circshift(bnd,[-1,0]); %Compute the boundary length of the structure under consideration
disp(['Boundary length = ' num2str(sum(sqrt(d(:,1).^2 + d(:,2).^2)),'%.1f')])

%% ----- DIPOLE PROPERTIES -----
r0 = [-4.5*R{1} + [0,rep(2)*R{2}(:,2)/2] + [0,.15],.125]; %position
polar = { [0,1,0] }; %polarization (y-polarized)

%% ----- AUXILIARY SQUARE MESH FOR PLOTTING -----
bndspc = .65; %boundary spacing

xstartend = [r0(1)-3.5,max(bnd(:,1)) + 1.25*bndspc];
ystartend = [max(bnd(:,2)) - 3.5498-sqrt(3)/2,max(bnd(:,2))+1.25*bndspc];
xyfrac = ((max(bnd(:,2))+1.25*bndspc) - (R{2}(:,2)*2) ) / (xstartend(2) - xstartend(1));
x = linspace(xstartend(1),xstartend(2),400);
y = linspace(ystartend(1),ystartend(2),ceil(300*xyfrac));

[X,Y] = meshgrid(x,y);
%% ----- INCLUSIONS -----
%inclusion basic shape
theta = linspace(0,2*pi,Ncirc); theta=theta(1:end-1).'+theta(2)/2;
circ = [cos(theta),sin(theta)]/(2*adivd);

%--- regions where we reduce the number of points around the circle ---
%entire bottom region
Ncircmod(1) = 20;
modregion{1}.x = [min(bnd(:,1)),min(bnd(:,1)),max(bnd(:,1)),max(bnd(:,1))]+[1,1,-1,-1]*.75;
modregion{1}.y = [bnd(1,2),bnd(2,2),bnd(2,2),bnd(1,2)]+[.75,-4.5,-4.5,.75];

%left, right and bottom edges of holes
Ncircmod(end+1:end+3) = 28;
modregion{end+1}.x = [bnd(1,1),bnd(1,1),max(bnd(:,1)),max(bnd(:,1))]; %bottom
modregion{end}.y =    bnd(1,2)+[0,.75,.75,0];
modregion{end+1}.x = bnd(1,1) + [0,0,.75,.75]; %left
modregion{end}.y = [bnd(1,2)+.75,max(bnd(:,2))-[.75,.75],bnd(1,2)+.75];
modregion{end+1}.x = max(bnd(:,1)) + [-.75,-.75,0,0]; %right
modregion{end}.y = [bnd(1,2)+.75,max(bnd(:,2))-[4.5,4.5],bnd(1,2)+.75];

%define a (local) hdata density for each modregion
hdatafac = [[1,1,1]*6,[1,1,1]*3]; %"multiplication factor" for these regions
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
    hdatamod{4}(x,y) + ...  %modified regions
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
    hold off; axis equal;
    xlim(minmax([minmax(X)+[-1,1]*bndspc/2,minmax(bnd(:,1))+[-1,1]*bndspc*1.25]))
    ylim(minmax([minmax(Y)+[-1,1]*bndspc/2,minmax(bnd(:,2))+[-1,1]*bndspc*1.25]))
    drawnow
end

%% ----- MESH THE ACTUAL STRUCTURE -----
[mesh.p,mesh.t,mesh.nodes,mesh.cnct] = geomFancyInclusion(bnd,incl,hdata,plotresults);

%% ----- EXTERNAL POTENTIAL DUE TO DIPOLE -----

phiext_spatial = zeros(size(mesh.p,1),numel(polar)); %preallocate
for pp = 1:numel(polar)
    phiext_spatial(:,pp) = calcDipoleSource([mesh.p,zeros(size(mesh.p,1),1)],r0,polar{pp});
end

if strcmpi(plotresults,'yes') || all(plotresults == 1)
    [cols,nams]=flatcolors(); hold on;
    plot(r0(:,1),r0(:,2),'p','MarkerSize',12.5,'MarkerFaceColor',cols{10}*.65+cols{6}*.35,'MarkerEdgeColor',cols{4}); hold on;
    drawnow;
end

%% ----- EXTERNAL POTENTIAL ON AUXILIARY SQUARE MESH -----

phiextsq_spatial = zeros(numel(X),numel(polar)); %preallocate
for pp = 1:numel(polar)
    phiextsq_spatial(:,pp) = calcDipoleSource([X(:),Y(:),zeros(numel(X),1)],r0,polar{pp});
end

if strcmpi(plotresults,'yes') || all(plotresults == 1)
    hold on;
    plot([min(x),max(x),max(x),min(x),min(x)],[min(y),min(y),max(y),max(y),min(y)],'-','Color',cols{6});
    drawnow;
    hold off;
end

%% ----- MODEL(S) SETUP -----
if ~exist('models','var') || isempty(models)
    ef_eV = .2; L = 400e-9;
    models = struct('type','isotropic', 'ef_eV',ef_eV,'gam_eV',6e-3,'L',L);
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
    models = struct('type','magneto','ef_eV',ef_eV,'L',L,'gam_eV',1e-3 ,'B',Bval);
    addstr = [addstr '_B' num2str(Bval)];
end

%Print model type
fprintf('%g setup-models are considered:\n',numel(models))
modelsprint = struct2table(models(:)); disp(modelsprint); clear modelsprint

%% ----- PRINT THE SAVE NAME MODIFIER IF IT EXISTS -------
if ~isempty(addstr)
    fprintf('Savename modifier is: %s\n\n',addstr);
end

%% ----- TEMPORAL ENVELOPE -----
%(NOTE: Right now, we do not allow more than one excitation "setting" -
% thus, all ef_eV in models must be IDENTICAL)

omegaf = unique(models.ef_eV)/hbar_eV;
fprintf('Temporal envelope setup:\n')
excitation.Omega = excitation.center_eV/(omegaf*hbar_eV); %Pulse energy-width (omegaF-normalized units)
fprintf('   Center frequency manually set to %.2f meV (= %.2f omegaF)\n',excitation.center_eV*1e3,excitation.Omega)
excitation.tw = 2*sqrt(log(2))/(excitation.width_eV/hbar_eV)*omegaf; %Pulse time-width (omegaF-normalized units)
excitation.tw_FWHM = excitation.tw*(4*sqrt(log(2))); %Pulse time-width FWHM (omegaF-normalized units)
fprintf('   Pulse-width manually set to tw = %.2f fs = %.2f 1/omegaF (~%.3f eV) \n',...
    excitation.tw/omegaf*1e15,excitation.tw,excitation.width_eV)

excitation.tdelay = 6*excitation.tw; %Is limited by a desire to let the FFT correspond to the true FT
%- should be large enough that field is very nearly zero for t = 0.
excitation.tdelay = (ceil(excitation.tdelay*excitation.Omega/pi-1/2)+1/2)*pi/excitation.Omega; %We correct the value of tdelay so
%that sin(Omega*tdelay) = 1, i.e. that the peak potential is unity (this is achieved by ensuring that tdelay = (n+1/2)pi/Omega -while still being larger than 8*tw)

if ~isfield(excitation,'profile') || strcmpi(excitation.profile,'pulse'); %Gaussian+sinusoidal
    phiext_temporal = @(T) sin(excitation.Omega*T).*exp(-(T-excitation.tdelay).^2/(4*excitation.tw^2));
    fprintf('      Sinusoidal limited by broad Gaussian (FWHM = %.1f fs, ~%.1f meV)\n\n',...
        excitation.tw_FWHM/omegaf*1e15,excitation.width_eV*1e3)
elseif strcmpi(excitation.profile,'cw'); %Let times t>tdelay be purely sinusoidal
    phiext_temporal = @(T) sin(excitation.Omega*T).*...
        (exp(-(T-excitation.tdelay).^2/(4*excitation.tw^2)).*double(T<=excitation.tdelay) ...
        + double(T>excitation.tdelay));
    fprintf('      CW-envelope chosen constant for t > t_delay \n')
    fprintf('      Sinusoidal limited by half-Gaussian (FWHM = %.1f fs, ~%.1f meV)\n\n',...
        excitation.tw_FWHM/omegaf*1e15,excitation.width_eV*1e3)
end

%Create the compound spatial * temporal phi-external function
phiext = @(T) bsxfun(@times,phiext_spatial,phiext_temporal(T).');
excitation.temporal_fstr = char(phiext_temporal); %save for possible postprocessing use

if ~isfield(excitation,'Tmax') %If Tmax not specified, we request 25 periods of Omega
    excitation.Tmax = [0, 3*excitation.tdelay + 10/excitation.Omega*2*pi, 2500]; %start-time, end-time, number of recorded steps (linearly separated)
end

if strcmpi(plotresults,'yes') || all(plotresults == 1) || strcmpi(plotresults,'excitation')
    set_figsize(5,10,10)
    Tplot = linspace(excitation.Tmax(1),excitation.Tmax(2),excitation.Tmax(3));
    plot(Tplot,phiext_temporal(Tplot),'-','color',[1,1,1]*.5); drawnow;
end


%% TIME-STEPPING
[T,rho,phiind] = timeEvolution(mesh,models,excitation.Tmax,phiext);
fprintf('%g time-steps, of mean size %.3g\n',numel(T),mean(diff(T)))
fprintf('phiind array is %g x %g, requiring %g GB of memory\n\n',size(phiind,1),size(phiind,2),getfield(whos('phiind'),'bytes')/1e9)

%% ----- SAVE DATA ON MESH -----
savevars = {'T','phiind','excitation','models','phiext_spatial','mesh',...
            'r0','polar','repincl','Ncirc','bndshape','bnd','incl'};
datapath = '../output/timestep/';
name = ['timestepDipole_' bndshape addstr];
save([datapath name],savevars{:},'-v7.3')
fprintf('Mesh-specific results (i.e. phiind) saved to %s\n\n',[datapath name '.mat'])

%% REFINE THE RESULTS AROUND THE ZOOM-IN MESH ON A SQUARE MESH

Vsq = calcCoulombGeneralPoints(mesh.p,mesh.t,[X(:),Y(:)]); 
fprintf('   Storing square-gridded Coulomb matrix Vsq requires %g GB memory\n',getfield(whos('Vsq'),'bytes')/1e9)

if numel(models) == 1; %phiind (and phiext) are implicitly normalized by the system size (models.L) -- regardless, the absolute values of phiindext are meaningless in the linear regime
    phiindsq = Vsq*rho/(4*pi*eps0);
else
    for mm = 1:numel(models)
        phiindsq{mm} = Vsq*rho{mm}/(4*pi*eps0);
    end
end
fprintf('phiindsq array is %g x %g, requiring %g GB of memory\n\n',size(phiindsq,1),size(phiindsq,2),getfield(whos('phiindsq'),'bytes')/1e9)
        
clear Vsq rho

%% ----- SAVE DATA ON SQUARE MESH-----
savevars = {'T','phiindsq','excitation','models','phiext_spatial','mesh',...
            'x','y','r0','polar','bndshape','bnd','incl'};
datapath = '../output/timestep/';
name = ['timestepDipole_' bndshape addstr '_SquareMesh'];
save([datapath name],savevars{:},'-v7.3')
fprintf('Square-mesh-specific results (i.e. phiindsq) saved to %s\n\n',[datapath name '.mat'])
