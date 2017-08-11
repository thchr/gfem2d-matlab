function runRibbonLatticeShortRange(latticetype,numcell,Ns,Ncirc,adivd,CoulombMethod,Nk,numeigs)

fprintf('\nCalculation (SHORT RANGE ONLY) commenced on | %s\n',datestr(now)); 
addpath(genpath('..'))

%%  ALLOW MULTITHREADING ON CLUSTER
fprintf('\n|----- MULTITHREADING SETUP -----|\n')
%SetMatlabMultithreading;

%% SETUP DEFAULTS
if nargin < 1 || isempty(latticetype)
    latticetype = 'triangular'; %Lattice type (square is the alternative)
end
if nargin < 2 || ~exist('numcell','var') || isempty(numcell)
    numcell = 2; %Number of repeated unit cells in ribbon (vertical dir.)
end
if nargin < 3 || ~exist('Ns','var') || isempty(Ns)
    Ns = 70; %Number of points along each unit cell side
end
if nargin < 4 || ~exist('Ncirc','var') || isempty(Ncirc) 
    Ncirc = 30; %Number of points used in circular inclusion
end 
if nargin < 5 || ~exist('adivd','var') || isempty(adivd) 
    adivd = 2; %Period divided by inclusion diameter
end  
if nargin < 6 || ~exist('CoulombMethod','var') || isempty(CoulombMethod)
    CoulombMethod = '2Dsupercell'; %Calculation method for Coulomb matrix (can be '1Dseries' or '2Dsupercell')
end 
if nargin < 7 || ~exist('Nk','var') || isempty(Nk)
    Nk = 25; %Number of k-points
end 
if nargin < 8 || ~exist('numeigs','var') || isempty(numeigs) 
    numeigs = 48; %Number of requested eigenmodes
end 

%Print setup
fprintf('\n|----- SETUP -----|\n')
strf = @(instr,strtype) [repmat(' ',1,14-numel(instr)) instr ' | %' strtype '\n']; %Convenient print string as function
fprintf([strf('lattice type','s'), strf('numcell','g'), strf('Ns','g'), strf('Ncirc','g'), strf('adivd','g'), strf('CoulombMethod','s'), strf('Nk','g'), strf('numeigs','g')],...
         latticetype             ,numcell             , Ns             , Ncirc           , adivd            , CoulombMethod            , Nk            , numeigs)

%% MESHING & MOMENTA
fprintf('\n|----- MESHING -----|\n')
theta = linspace(0,2*pi,Ncirc+1); theta = theta(1:end-1); %Inclusion boundary
inclusion = 1/(2*adivd)*[cos(theta);sin(theta)].';        %----||----
blochmesh = geomRibbonInclusion(latticetype, [1,numcell],inclusion,Ns,@(x,y) .015*250/Ns*ones(size(x)),0);

%Momentum list
kvec = bsxfun(@times,linspace(0,1/2,Nk)',blochmesh.Grib); 
maxk = max(sqrt(kvec(:,1).^2 + kvec(:,2).^2));

erfccut = 5.5; %Cutoff; such that erfc(erfccut) has a value less than 1e-14
%% ACTUAL CALCULATION AND PLOTTING
zeta = zeros(numeigs,Nk); eigrho = zeros(size(blochmesh.p,1),numeigs,Nk); eigphi = eigrho; %Preallocate

fprintf('|----- DISPERSION CALCULATION (%g k-points) -----|\n',Nk)
tick = tic;
for kk = Nk:-1:1
    %--- ACTUAL CALCULATION ---
    switch CoulombMethod
        case '1Dseries'
            V = calcRibbonCoulombViaSeries(blochmesh,kvec(kk,:),10);
        case '2Dsupercell'
            if kk == Nk; 
                [~,A_uc] = calcReciprocal(blochmesh.R); E = sqrt(pi/A_uc);
                blochmesh.R{2} = 3*numcell*blochmesh.R{2}; 
                blochmesh.G = calcReciprocal(blochmesh.R);  
            end
            nmax(1,1) = max(ceil(erfccut/(E*norm(blochmesh.R{1},2))+1),2);
            nmax(1,2) = max(ceil(erfccut/(E*norm(blochmesh.R{2},2))+1),2);
            
            V = calcLatticeCoulombShortRange(blochmesh,kvec(kk,:),nmax,E);
    end
    [DRi,DL] = calcLatticeDifferentialIsotropic(blochmesh,kvec(kk,:));
    [zeta(:,kk),eigrho(:,:,kk)] = calcEigen(V,DL,DRi,numeigs,[],true);
    eigphi(:,:,kk) = V*eigrho(:,:,kk);
    
    %--- PRINT PROGRESS ---
    fprintf('--- DISPERSION LOOP %g/%g (%.2f/%.2f min) ---\n\n',Nk+1-kk,Nk,toc(tick)/60,toc(tick)/60*Nk/(Nk+1-kk));
end

%% SAVE DATA
clear W DRi DL %Don't want to save this; it's heavy
savedir = '../output/ribbons/shortrangeonly/';
savename = [latticetype '_numcell' num2str(numcell) '_' CoulombMethod '_Ns' num2str(Ns) '_Ncirc' num2str(Ncirc)];
save([savedir savename '.mat'])

%Print final messages
fprintf('\n\n All data saved to   | %s\n\n', [savedir savename])
fprintf('Calculated finished on   | %s\n  ', datestr(now))