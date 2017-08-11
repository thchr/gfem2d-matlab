function runDispersionAntidotLattice(adivd,Ne,Ncirc,Nk,latticetype,meshfun,plotgraphs)

%Load necessary subfolders and constants
addpath(genpath('..'));
ConstantsUnits0; 

%Check input
if ~exist('adivd','var') || isempty(adivd)
	adivd = 2; %Must always be larger than 1! Radius is L/(2*adivd) where a = L
end
if ~exist('Ne','var') || isempty(Ne)
	Ne = 120;
end
if ~exist('Ncirc','var') || isempty(Ncirc)
	Ncirc = 150;
end
if ~exist('Nk','var') || isempty(Nk)
	Nk = 35;
end
if ~exist('latticetype','var') || isempty(latticetype)
	latticetype = 'triangular';
end
if ~exist('meshfun','var') || isempty(meshfun)
	meshfun = 'antidot';
end
if ~exist('plotgraphs','var') || isempty(plotgraphs)
    plotgraphs = 0;
end

%Inclusion boundary
theta = linspace(0,2*pi,Ncirc+1); theta = theta(1:end-1);
inclusion = 1/(2*adivd)*[cos(theta);sin(theta)].';

%Meshing
blochmesh = geomPeriodicInclusion(latticetype,inclusion,Ne,meshfun,plotgraphs);

%Momentum list
[kvec,kplot,kmark] = irreducibleFBZ(latticetype,Nk,plotgraphs);

%parpool;
ticEL=tic;
for kk = 1:size(kvec,1);
    zeta(:,kk) = calcLatticeEigen(blochmesh,kvec(kk,:));
    fprintf('\n\n     DISPERSION LOOP: %g/%g (%.1f min/%.1f min)\n\n\n',kk,size(kvec,1),toc(ticEL)/60,size(kvec,1)/kk*toc(ticEL)/60);
end

%% SAVE DATA
savename = ['Dispersion' upper(latticetype(1)) latticetype(2:end) ...
            'AntidotLattice_adivd' num2str(adivd) ...
            '_Ne' num2str(Ne) 'Nc' num2str(Ncirc) ...
            '_Meshfun' upper(meshfun(1)) meshfun(2:end)];

%Save data
save(['../output/' savename '.mat']);

%% PLOT SQUARE ROOT OF ZETA VALUES (EIGENVALUES)

if exist('plotgraphs','var') && plotgraphs == 1
	R = blochmesh.R;
	G = calcReciprocal(R);
	
	[cols,nams]=flatcolors;
	set_figsize([],14,14); hold on;

	for mm = 2:size(kmark.n,1)-1
		plot(kplot(kmark.n(mm,1))*[1,1],[0,1]*7,'-','Color',[1,1,1]*.25,'LineWidth',.45)
	end
	
	plot(kplot, sqrt(zeta),'-','color',cols{2},'LineWidth',1.5)
	for n = -5:5
		for m = -5:5
			Gnm = G{1}*n + G{2}*m;
			
			val = sqrt( (kvec(:,1) + Gnm(1)).^2 + (kvec(:,2) + Gnm(2)).^2);
			
			%drawnow
			plot(kplot,sqrt(val),':','color',cols{4});
		end
	end
	hold off
	box on
	set(gca,'Layer','Top','Fontsize',8)
	format_ticks(gca,{kmark.symbol{kmark.n(:,2)}},[],[kplot(kmark.n(:,1))],[],[],[],0.005)
	axis tight
	ylim([0,3])
	xlim([0,max(kplot)])
	xlabh = get(gca,'XLabel');
	set(xlabh,'Position',get(xlabh,'Position') - [0 .1 0])
	xlabel('Momentum','Fontsize',10)
	ylabel('Dimensionless resonance value','Fontsize',10)

    %Export figure if opened
	export_fig(['../figures/' savename],'-pdf')
end