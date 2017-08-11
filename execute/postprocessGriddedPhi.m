function postprocessGriddedPhi(kchoice,Nx)

fprintf('\nCalculation commenced on | %s\n',datestr(now)); 
addpath(genpath('..'))

%%  ALLOW MULTITHREADING ON CLUSTER
fprintf('\n|----- MULTITHREADING SETUP -----|\n')
SetMatlabMultithreading;

%% PREALLOCATION
if ~exist('kchoice','var') || isempty(kchoice)
    kchoice = 12;
end
if ~exist('Nx','var') || isempty(Nx)
    Nx = 40;
end

%% LOAD MODE PROFILE DATA
load(['../output/ribbons/scissored_eigvec/scmag_triangular_numcell16_cut40_Ns53_Ncirc36_savephi1_kpnt' num2str(kchoice) 'd18.mat']);
remesh = mesh.remesh;
%Some functions are somehow in trouble, and cannot be called, despite being
%defined; we hardcode-redefine them here:
remesh.values = @(f)[f;f(remesh.partner(:,1))];

%% MAKE A NEW QUADRATIC MESH
ymax = max(remesh.p(:,2));
xlims = .25+3.5+[-.5,+.5];
ylims  = ymax + [-4,1];

Ny = round(Nx*diff(ylims)/diff(xlims));

x = linspace(xlims(1),xlims(2),Nx+1); x = x(1:end-1) + diff(x(1:2))/2;
y = linspace(ylims(1),ylims(2),Ny+1); y = y(1:end-1) + diff(y(1:2))/2;
[X,Y] = meshgrid(x,y);

fprintf('--- A SQUARE MESH WAS CREATED (Nx = %g, Ny = %g, Ntotal = %g) ---\n',Nx,Ny,Nx*Ny)

%% CONSTRUCT AUXILIARY MESH WHICH LIES IN THE EDGE UNIT CELL

M = 25;  %Factor with which the BZ is multiplied (should be greater than 8 for a 16-ribbon)
extendbz.x = [-M,-M,M,M]*mesh.R{2}(1) + [-1,1,1,-1]*mesh.R{1}(1)/2;
extendbz.y = [-M,-M,M,M]*mesh.R{2}(2) + [-1,1,1,-1]*mesh.R{1}(2)/2;

fprintf('--- CONSTRUCTING AUXILLIARY MESH ---\n')
maxshift = 10; %Hardcoded top-limit on shifting back and forth
aX = X(:); aY = Y(:); loopbreak = 1;
for nn = 1:numel(aX);
    while loopbreak < 2*maxshift;
        if inpolygon(aX(nn),aY(nn),extendbz.x,extendbz.y)
            loopbreak = 1;
            break
        else
            if loopbreak <= maxshift;
                aX(nn) = aX(nn)-mesh.Rrib(1);
            elseif loopbreak == maxshift+1;
                aX(nn) = X(nn)+mesh.Rrib(1);
            elseif loopbreak > maxshift+1
                aX(nn) = aX(nn)+mesh.Rrib(1);
            end
            loopbreak = loopbreak + 1;
        end
    end
    if mod(nn,1000) == 0;
        fprintf('   %g/%g\n',nn,numel(aX))
    end
end
ar = {aX,aY};

%% COMPUTE NEW POTENTIALS AND ASSIGN THEM TO CELL ARRAY
numeigs = 100;
phixy = cell(3,2,numeigs);
for ksgn = [-1,+1];
    ksgni = (ksgn == -1)*1 + (ksgn == +1)*2;
    fprintf('\n\n--- CALCULATING GRIDDED COULOMB OPERATOR (sign of kx = %g) ---\n',ksgn)
    %Calculate the Coulomb operator for the new grid
    V = calcRibbonCoulombViaSeries(mesh,kvec(ksgni,:),bloch.n_genexpint,ar);
    
    %Calculate the new gridded potentials
    for mm = 1:3;
        phixy_matrix = V*rhoeig{mm}{ksgni};
        for bb = 1:numeigs;
            phixy{mm}{ksgni}{bb} = reshape(phixy_matrix(:,bb),size(X));
        end
    end
end

%% SAVE DATA
savevars = {'x','y','Nx','Ny','X','Y','kvec','bloch','xlims','ylims','kchoice',...
            'phixy','mesh','Ns','Ncirc','numcell','cut','adivd','latticetype',...
            'models','eneeig_eV','ef_eV','L'};

savedir = '../output/ribbons/scissored_eigvec/';
savename = ['phi_gridNx' num2str(Nx) '_numcell16_cut40_Ns53_Ncirc36_kpnt' num2str(kchoice) 'd18.mat'];
savepath = [savedir savename]; 

save(savepath,savevars{:})

%Print final messages
fprintf('\n\n All data saved to   | %s\n\n', savepath)
fprintf('Calculated finished on   | %s\n  ', datestr(now))
%%
% PLOT NEW XY-GRIDDED POTENTIAL (OLD NOT WORKING)
% if 0
%     set_figsize(2,20,25)
%     subaxis(1,2,1)
%     imagesc(x,y,real(phixy)); set(gca,'YDir','normal')
%     caxis([-1,1]*max(abs(real(phixy(:)))));
%     colormap(bluewhitered);
%     axis equal
%     xlim(minmax(x)); ylim(minmax(y));
%     
%     subaxis(1,2,2)
%     phimesh = real(remesh.values(phieig{MM}{ksgni}(:,BB)));
%     for rr = -5:5
%         trisurf(remesh.t,remesh.p(:,1)+rr*mesh.Rrib(1),remesh.p(:,2),phimesh,'EdgeColor','none'); hold on
%     end
%     view(2);
%     caxis([-1,1]*max(abs(phimesh)));
%     colormap(bluewhitered);
%     axis equal
%     xlim(minmax(x)); ylim(minmax(y));
% end