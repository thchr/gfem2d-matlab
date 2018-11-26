function runWilsonLoop(latticetype,Ns,Ncirc,adivd,loopkN,looptype,ind,addstr,plotfigs,savefig)

fprintf('\n2D lattice calculation commenced on | %s\n',datestr(now));
addpath(genpath('..'))

%%  ALLOW MULTITHREADING ON CLUSTER
fprintf('\n|----- MULTITHREADING SETUP -----|\n'); nopool = 1; %Don't open a multiworker pool
SetMatlabMultithreading;

%% SETUP DEFAULTS
if nargin < 1 || isempty(latticetype)
    latticetype = 'triangular'; %Lattice type (square is the alternative)
end
if nargin < 2 || ~exist('Ns','var') || isempty(Ns)
    Ns = 70; %Number of points along each unit cell side
end
if nargin < 3 || ~exist('Ncirc','var') || isempty(Ncirc)
    Ncirc = 36; %Number of points used in circular inclusion
end
if nargin < 4 || ~exist('adivd','var') || isempty(adivd)
    adivd = 2; %Period divided by inclusion diameter
end
if nargin < 5 || ~exist('loopkN','var') || isempty(loopkN)
    loopkN = 30; %# of k-points in noncontractible loop
end
if nargin < 6 || ~exist('looptype','var') || isempty(looptype)
    looptype = 'xA'; %select a predefined loop type (xA & xB implemented now; yA & yB could also be implemented)
end
if nargin < 7 || ~exist('ind','var') || isempty(ind)
    ind = {[1,2,3],[4,5]}; %Multiplet/multi-band grouping of Berry phases
end
if nargin < 8 || ~exist('addstr','var')
    addstr = []; %Whether to add a modifier string to the savename ([] = nothing added)
elseif ~isempty(addstr) && addstr(1)~='_'
    addstr = ['_' addstr];
end
if nargin < 9 || ~exist('plotfigs','var') || isempty(plotfigs)
    plotfigs = 0; %Whether to plot a figure associated with the Berry phases or not
end
if nargin < 10 || ~exist('savefig','var') || isempty(savefig)
    savefig = 0; %Whether to save this figure
end

%% PRINT SETUP

strf = @(instr,strtype) [repmat(' ',1,14-numel(instr)) instr ' | %' strtype '\n']; %Convenient print string as function

fprintf('\n|----- SETUP -----|\n')
fprintf([strf('lattice type','s'), strf('Ns','g'), strf('Ncirc','g'), strf('adivd','g'), strf('loopkN','g'), strf('looptype','s'), strf('ind','s'),                                                         strf('addstr','s')],...
         latticetype             , Ns            , Ncirc            , adivd            , loopkN            , looptype,             strjoin(cellfun(@(x) ['[' num2str(x) '], '],ind,'UniformOutput',false)), addstr)

%% MESHING
fprintf('\n|----- MESHING UNIT CELL -----|\n')
theta = linspace(0,2*pi,Ncirc+1); theta = theta(1:end-1); %Inclusion boundary
inclusion = 1/(2*adivd)*[cos(theta);sin(theta)].';        %----||----
mesh = geomPeriodicInclusion(latticetype,inclusion,Ns,@(x,y) .015*250/Ns*ones(size(x)),plotfigs);

%% MOMENTA
G = calcReciprocal(mesh.R);
Nk = 31;
switch looptype
    case 'xA' %An 'up-down' type loop along ky, with each loop having fixed kx
        kstart = [linspace(-G{1}(1),G{1}(1),Nk) ;...
                  linspace(-G{2}(2),-G{2}(2),Nk)].'/2 ;
        loop.G = G{2};
    case 'yA' %A 'left-right' type loop along kx, with each loop having fixed ky
        kstart = [linspace(-G{1}(1),-G{1}(1),Nk) ;
                  linspace(-G{1}(2), G{1}(2), Nk)/2 ].';
        loop.G = 2*G{1} + G{2};
    otherwise
        error('looptype must be one of a set of predefined looptype-strings: ''xA'' or ''yA''')           
end

loop.N = loopkN;

%% WILSON LOOP
%prepare for plotting, if requested 
if plotfigs==1
    cols=flatcolors;
    figberry=set_figsize([],10,17); figloop=set_figsize([],15,15); 
    [~,~,~,BZ] = irreducibleFBZ(latticetype,5,0); 
    patch([BZ(:,1);BZ(1,1)],[BZ(:,2);BZ(1,2)],cols{8}*.15+.85,...
                              'LineStyle','-','EdgeColor',cols{4}*.85+.15);
else 
    figloop = 0; 
end
%preallocate results-vector
berry = cell(numel(ind),1);
for ii = 1:numel(ind)
    berry{ii} = zeros(numel(ind{ii}),size(kstart,1));
end
wc = berry;

fprintf('\n\n|----- BERRY PHASES OF NONCONTRACTIBLE LOOPS (%g distinct starting k-points) -----|\n',size(kstart,1))
tkpoints=tic;
for kkstart = 1:size(kstart,1)
    fprintf('\nStarting k-point #: %g/%g\n',kkstart,size(kstart,1))
    loop.kstart = kstart(kkstart,:);
    
    %%
    
    [berrytemp, wctemp] = WilsonLoop(loop,mesh,ind,figloop);
    for ii = 1:numel(ind)
        berry{ii}(:,kkstart) = berrytemp{ii};
        wc{ii}(:,kkstart) = wctemp{ii};
    end
    if plotfigs == 1
        figure(figberry);
        plotkindex = int32(looptype(1))-119; %the above is a bit of a hacky solution to choose 1st or 2nd column if x (=120) or y (=121) are first characters of plottype
        for ii = 1:numel(ind)
            subaxis(numel(ind),1,numel(ind)+1-ii,'ML',.165,'MB',.085,'MT',.04,'MR',.025,'S',0.035)
            
            plot([0,0],[-1,1]*pi,'-','color',cols{8}); hold on
            plot(kstart(1:kkstart,plotkindex),berry{ii}(:,1:kkstart),'o-k','markerfacecolor',cols{ii},'color',cols{ii}); hold off;
            
            if ii == 1
                xlabel('Wavevector \itk_xa')
                set(gca,'XTick',[-1,0,1]*pi,'XTickLabel',{'-\pi','0','\pi'})
            else
                set(gca,'XTick',[-1,0,1]*pi,'XTickLabel',[])
                title(['\ita/d = ' num2str(adivd)],'FontWeight','Normal')
            end
        
            xlim(minmax(kstart(:,plotkindex))); ylim([-1,1]*pi);
            ylabel('Berry phase \it\phi_n');
            set(gca,'YTick',[-1:.5:1]*pi,'YTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'},...
                'TitleFontWeight','normal','TickDir','out','TickLength',[1,1]*0.015); %'TitleFontSizeMultiplier',1,
            text(pi*.925,pi,['Bands ' listEnglishGrammar(ind{ii})],'Margin',1,...
                'HorizontalAlignment','Right','BackgroundColor','w','Color',cols{ii})
            
            hold off
        end
        drawnow;
        
        if savefig == 1
            export_fig(['../figures/qshe/BerryPhases_' latticetype '_adivd' num2str(adivd) '_Ns' num2str(Ns) '_Ncirc' num2str(Ncirc)],'-pdf')
        end
    end
    fprintf('   Total time spent: %.1f min (/~%.1f min)\n',toc(tkpoints)/60,size(kstart,1)/kkstart*toc(tkpoints)/60)
end
savepath = ['../output/qshe/BerryPhases_' latticetype '_adivd' num2str(adivd,'%.1f') ...
            '_Ns' num2str(Ns) '_Ncirc' num2str(Ncirc) '_Nkloop' num2str(loopkN) '_' looptype addstr '.mat'];
save(savepath)

fprintf('\n\n|----- ALL DATA SAVED TO: %s -----|\n\n',savepath)
