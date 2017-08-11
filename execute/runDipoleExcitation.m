function runDipoleExcitation(bndshape,auxmesh,addstr,models,repincl,R,adivd,cut,Ncirc,plotresults,usedefect)

addpath(genpath('../'));
ConstantsUnits0;

%%  ALLOW MULTITHREADING ON CLUSTER
fprintf('\n|----- MULTITHREADING SETUP -----|\n'); nopool = 1; %Don't open a multiworker pool
SetMatlabMultithreading;

%% DEFAULT INPUT
if ~exist('bndshape','var') || isempty(bndshape);
    bndshape = 'rectanglewedge'; %See posibilities under the switch statement for bndshape
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
if ~exist('usedefect','var') || isempty(usedefect); %Basically an obsolete functionality
    usedefect = 'no';
    collisiontest = 0;
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
[p,t,nodes,cnct] = geomFancyInclusion(bnd,incl,[],plotresults);

%% ----- EXTERNAL POTENTIAL -----
r0 = [-4.5*R{1} + [0,rep(2)*R{2}(:,2)/2] + [0,.15],.125]; %position
%polar = {[1,0,0],[0,1,0],[0,0,1]};
polar = { [0,1,0] }; 
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
        xstartend = [r0(1)-3.85*bndspc,max(p(:,1))-4.075*bndspc];
        xyfrac = ((max(p(:,2))+1.25*bndspc) - (R{2}(:,2)*2) ) / (xstartend(2) - xstartend(1));
        x = linspace(xstartend(1),xstartend(2),400); y = linspace(R{2}(:,2)*2,max(p(:,2)+1.25*bndspc),ceil(300*xyfrac));
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
    ef_eV = 0.2; L = 400e-9;
    models(1,1) = struct('ef_eV',ef_eV,'L',L,'gam_eV',1e-3 ,'ene_eV',60.3e-3,'B',2);
    models(1,2) = struct('ef_eV',ef_eV,'L',L,'gam_eV',.5e-3,'ene_eV',60.3e-3,'B',2);
    models(2,1) = struct('ef_eV',ef_eV,'L',L,'gam_eV',1e-3 ,'ene_eV',60.3e-3,'B',4);
    models(2,2) = struct('ef_eV',ef_eV,'L',L,'gam_eV',.5e-3,'ene_eV',60.3e-3,'B',4);
    models(3,1) = struct('ef_eV',ef_eV,'L',L,'gam_eV',1e-3 ,'ene_eV',60.3e-3,'B',8);
    models(3,2) = struct('ef_eV',ef_eV,'L',L,'gam_eV',.5e-3,'ene_eV',60.3e-3,'B',8);
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
    models(1) = struct('ef_eV',ef_eV,'L',L,'gam_eV',1e-3 ,'ene_eV',60.3e-3,'B',Bval);
    models(2) = struct('ef_eV',ef_eV,'L',L,'gam_eV',.5e-3,'ene_eV',60.3e-3,'B',Bval);
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
V = calcCoulomb(p,t); whosV = whos('V');
fprintf('   Storing Coulomb matrix V requires %g GB memory\n\n',whosV.bytes/1e9)
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
    
    ticslv = tic; %We write it in this awful lengthy way to not eat up more memory than necessary
    results(mm).rho = ( (DL+1i*results(mm).sigma0/(4*pi*eps0*models(mm).L*omega)*DR*V) ) \ ...
                      ( -1i*results(mm).sigma0/(omega*models(mm).L^2)*DR*phiext );
    fprintf('   Matrix system solve in %.2f min\n\n',toc(ticslv)/60)
    
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
name = ['dipoleExcite_' bndshape '_' auxmesh addstr];
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