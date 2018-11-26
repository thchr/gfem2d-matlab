clear; clc; close all;
cols=flatcolors();
addpath(genpath('..'));
ConstantsUnits0;
%% Make a disk and calculate the polarizability and eigenstruct values
geom = 'disk';
switch geom
    case 'disk'
        hdata.fun = @(x,y) real(.0375*(cos(pi/2*(x.^2+y.^2)).^.5+.45).^(.75));
        [mesh.p,mesh.t,mesh.edge,mesh.edgecnct] = geomDisk(451,hdata,0); %Mesh
        sympairs = {[1,2],[12,13],[27,28],[48,49],[71,72],[102,103]};  %Choices of eigenindices corresponding to first 6 dipolar modes
        %Last pair should be [100,101] for really great meshes and [102,103] otherwise
    case 'square'
        [mesh.p,mesh.t,mesh.edge,mesh.edgecnct] = geomSquare(1,125,.021,0);
        sympairs = {[1,2],[6,7],[9,10],[14,15],[21,22],[32,33],[38,39]};
    case 'triangle'
        [mesh.p,mesh.t,mesh.edge,mesh.edgecnct] = geomEquiTriangle(1,151,.015,0); %Set side to 151 for a good mesh
        sympairs = {[1,2],[4,5],[6,7],[10,11],[13,14],[17,18]};
    case 'ellipse'
        mesh = geomEllipse([1,.5],501,.025,0);
        sympairs = [];
end
if isfield(mesh,'areavec'); error('areavec initialized too soon; danger'); end
mesh.edge = mesh.edge/sqrt(sum(meshArea(mesh.p,mesh.t)));  %Normalize to unit area
mesh.p = mesh.p/sqrt(sum(meshArea(mesh.p,mesh.t)));        %--||--

plotWireMesh(mesh.p,mesh.t,1)


%% Calculate polarizability and eigen-properties
[alpha,eigstruct] = calcPolarizability(mesh,110,1,sympairs); %Actual calculation

%%
% Print the values of modes with substantial oscillator strength
toldisp = 1e-3; nplot = [];
fprintf('   n  |     zeta    osc.x    osc.y\n');
for n = 1:size(eigstruct.zeta)
    if any(eigstruct.oscstrength(n,:)>toldisp)
        fprintf(' %3g  |  ',n)
        fprintf('%7.4f   ',eigstruct.zeta(n)); %Eigenvalue (forced equal pairwise in calcPolarizability)
        fprintf('%6.4f   %6.4f\n',eigstruct.oscstrength(n,:)); %Sum the values in each pair to get x- and y- oscillator strengths
        nplot(end+1) = n;
    end
end
fprintf('\n\n')

% Print the values of the first 6 dipolar modes
if ~isempty(sympairs)
    n = 1;
    fprintf(' n  |     zeta    osc.x    osc.y\n');
    for pair = sympairs
        fprintf(' %g  |  ',n)
        fprintf('%7.4f   ',eigstruct.zeta(pair{:}(1))); %Eigenvalue (forced equal pairwise in calcPolarizability)
        fprintf('%6.4f   %6.4f\n',sum(eigstruct.oscstrength(pair{:},:),1)); %Sum the values in each pair to get x- and y- oscillator strengths
        n = n+1;
    end
    nplot = [sympairs{:}];
end

%% New Coulomb matrix
Nmean = 225; delim = 0.5;
Nx = round(diff(minmax(mesh.p(:,1))+[-1,1]*delim)./diff(minmax(mesh.p(:,2))+[-1,1]*delim)*Nmean);
Ny = round(diff(minmax(mesh.p(:,2))+[-1,1]*delim)./diff(minmax(mesh.p(:,1))+[-1,1]*delim)*Nmean);
x = linspace(min(mesh.p(:,1))-delim,max(mesh.p(:,1))+delim,Nx);
y = linspace(min(mesh.p(:,2))-delim,max(mesh.p(:,2))+delim,Ny);
[X,Y] = meshgrid(x,y);

Vsq = calcCoulombGeneralPoints(mesh.p,mesh.t,[X(:),Y(:)]);

%% Plot associated modal profiles
set_figsize(2,40,25.5)
if isempty(sympairs)
    % Plot without symmetry constraints

    for nn = 1:min(numel(nplot),12)
        osc = eigstruct.oscstrength(nplot(nn),:);
        subaxis(3,4,nn,'M',.01,'S',.05,'MT',.05)
        %trisurf(mesh.t,mesh.p(:,1),mesh.p(:,2),eigstruct.phi(:,nplot(nn)),'EdgeColor','none'); shading interp;
        contourf(X,Y,reshape(Vsq*eigstruct.rho(:,nplot(nn)),Ny,Nx),16,'LineWidth',.5,'LineColor',cols{8}*.4+cols{4}*.6);
        hold on
        plot(mesh.edge(mesh.edgecnct(1:2:end),1),mesh.edge(mesh.edgecnct(1:2:end),2),'-k')
        hold off
       
        colormap(gca,'bluewhitered_mod')
        axis equal off;
        title(['\itn\rm = ' num2str(nplot(nn)) ...
               ', \zeta = ' num2str(eigstruct.zeta(nplot(nn)),'%.4f') ...
               ', \Delta\it_x\rm = ' num2str(round(osc(1)*1e4)/1e4,'%.4g') ...
               ', \Delta\it_y\rm = ' num2str(round(osc(2)*1e4)/1e4,'%.4g')], ...
               'Fontsize',10,'FontWeight','normal')
        drawnow; 
    end
else
    
    % Plot with symmetry constraints

    for nn = 1:min(numel(sympairs),6)
        pair = sympairs(nn);
        
        pairrho = eigstruct.rho(:,pair{:});
        pairphi = eigstruct.phi(:,pair{:});
        pairzeta = mean(eigstruct.zeta(pair{:}));
        for mm = 1:2
            ampl(1,mm) = integrateMeshFunction(mesh.p,mesh.t,mesh.p(:,1).*pairrho(:,mm));
            ampl(2,mm) = integrateMeshFunction(mesh.p,mesh.t,mesh.p(:,2).*pairrho(:,mm));
        end
        for ii = 1:2 %Plot x- and y-polarized mode separately (combine modes to get this polarization)
            rho = zeros(size(pairrho(:,1))); phi = zeros(size(pairphi(:,1)));
            for mm = 1:2
                rho = rho + pairrho(:,mm)*ampl(ii,mm);
                phi = phi + pairphi(:,mm)*ampl(ii,mm);
            end
            %Normalize eigenvectors such that <rho_n|V|rho_n> = 2*pi*zeta_n (consistent
            %with choice in your thesis, in our present notation)
            overlap = integrateMeshFunction(mesh.p,mesh.t,rho.*phi);
            rho = bsxfun(@times,rho,sqrt(2*pi*pairzeta./overlap).');
            %Perform analytical integration over each triangular element to obtain
            %normalized oscillator strengths along x and y (stored in 1st and 2nd column, respectively)
            osc(1) = abs(integrateMeshFunction(mesh.p,mesh.t,mesh.p(:,1).*rho)).^2;
            osc(2) = abs(integrateMeshFunction(mesh.p,mesh.t,mesh.p(:,2).*rho)).^2;
            
            phisq = Vsq*rho;
            subaxis(3,4,2*(nn-1)+ii,'M',.01,'S',.05,'MT',.05)
            contourf(X,Y,reshape(phisq,Ny,Nx),16,'LineWidth',.5,'LineColor',cols{8}*.4+cols{4}*.6);
            hold on
            plot(mesh.edge(mesh.edgecnct([1:end,1],1),1),mesh.edge(mesh.edgecnct([1:end,1],1),2),'-k','linewidth',1.25)
            hold off
          
            colormap(gca,'bluewhitered_mod')
            axis equal off;
            title(['\itn\rm = ' num2str(pair{:}(ii)) ...
                ', \zeta = ' num2str(pairzeta,'%.3f') ...
                ', \Delta\it_x\rm = ' num2str(round(osc(1)*1e4)/1e4,'%.4g') ...
                ', \Delta\it_y\rm = ' num2str(round(osc(2)*1e4)/1e4,'%.4g')], ...
                'Fontsize',11,'FontWeight','normal')
            drawnow; 
            
        end
    end
end
export_fig(['../figures/eigenproperties/Eigenproperties_' geom '_unitarea'],'-pdf','-opengl','-r400')
%% Plot the absorption cross-section for a set of material values
L=50e-9; ef_eV = 0.5; gam_eV = 1e-3;
omega = linspace(.05,.55,5000)/hbar_eV;

sigma = conducLRA(omega*hbar_eV,ef_eV,gam_eV);
zetaomega = 2i*eps0*omega*L./sigma;
Q_abs.x = imag(alpha.x(zetaomega*sqrt(pi),sqrt(pi)*L)).*omega/c/(pi*L^2); %Note the factors of sqrt(pi) here; this is because we are calculating 
Q_abs.y = imag(alpha.y(zetaomega*sqrt(pi),sqrt(pi)*L)).*omega/c/(pi*L^2); %for a structure of L = sqrt(area) = unity rather than L = radius = unity.


figure; set(gcf,'color','w') %Plot results
semilogy(omega*hbar_eV,Q_abs.x,'-','color',[.1,.1,.8],'linewidth',1.25); hold on; %x-cross-sectional efficiency
semilogy(omega*hbar_eV,Q_abs.y,'--','color',[.8,.1,.1],'linewidth',1.25); hold on; %y-cross-sectional efficiency
semilogy(omega*hbar_eV,omega/c.*imag( 2*L^3 * ( 2.8912 ./ (1.0977 - zetaomega) + ... %Analytical comparison
    0.1120 ./ (4.9140 - zetaomega) + ...
    0.0424 ./ (8.1337 - zetaomega) + ...
    0.0224 ./ (11.3079 - zetaomega) + ...
    0.0140 ./ (14.4675 - zetaomega) + ...
    0.0096 ./ (17.6205 - zetaomega) ) )/(pi*L^2),'-.k');

%misc clean-up of plot
xlim(minmax(omega*hbar_eV)); ylim(minmax(Q_abs.x)+[-1e-5,4.5])
legend({'Numerical; \itx\rm-polarization','Numerical: \ity\rm-pol.','Analytical: Table 3'},'Fontsize',9); legend boxoff
xlabel('Frequency (eV)','Fontsize',13); ylabel('Absorption efficiency','Fontsize',13);
title('\rmDisk (\itE\rm_F = 0.5 eV, \itR\rm = 50 nm, #vertices ~ 2300)','Fontsize',13)
set(gca,'Fontsize',8,'LineWidth',.5,'TickDir','out')
%export_fig('DiskAbsorptionEfficiency','-pdf')

