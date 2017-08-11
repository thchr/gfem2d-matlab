clear; clc; close all; 

addpath(genpath('..'));
ConstantsUnits0;
%% Make a disk and calculate the polarizability and eigenstruct values
[mesh.p,mesh.t] = geomDisk(1500,[],1); %Mesh
%[mesh.p,mesh.t] = geomDisk(70,.125,1); %Mesh

sympairs = {[1,2],[12,13],[27,28],[48,49],[71,72],[102,103]};  %Choices of eigenindices corresponding to first 6 dipolar modes
                                                               %Last pair should be [100,101] for really great meshes and [102,103] otherwise
[alpha,eigstruct] = calcPolarizability(mesh,110,1,sympairs); %Actual calculation

% Print the values of the first 6 dipolar modes
n = 1; 
fprintf(' n  |     zeta    osc.x    osc.y\n');
for pair = sympairs
fprintf(' %g  |  ',n)
fprintf('%7.4f   ',eigstruct.zeta(pair{:}(1))); %Eigenvalue (forced equal pairwise in calcPolarizability)
fprintf('%6.4f   %6.4f\n',sum(eigstruct.oscstrength(pair{:},:),1)); %Sum the values in each pair to get x- and y- oscillator strengths
n = n+1; 
end

%% Plot the absorption cross-section for a set of material values
L=50e-9; ef_eV = 0.5; gam_eV = 1e-3;
omega = linspace(.05,.55,5000)/hbar_eV;

sigma = conducLRA(omega*hbar_eV,ef_eV,gam_eV);
zetaomega = 2i*eps0*omega*L./sigma;
Q_abs.x = imag(alpha.x(zetaomega,L)).*omega/c/(pi*L^2);
Q_abs.y = imag(alpha.y(zetaomega,L)).*omega/c/(pi*L^2);


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

%% Misc stuff to test the overlap between various states
N = 15; %Number of states to test (from smallest zetaeig and upward)

%normalize self-overlap to one for easy comparison
diagoverlap = integrateMeshFunction(mesh.p,mesh.t,eigstruct.rho(:,1:N).*eigstruct.phi(:,1:N));
eigstruct.rho(:,1:N) = bsxfun(@times,eigstruct.rho(:,1:N),sqrt(1./diagoverlap).'); 
eigstruct.phi(:,1:N) = bsxfun(@times,eigstruct.phi(:,1:N),sqrt(1./diagoverlap).');

%calculate mutual overlap
overlap=zeros(N,N);
for ii = 1:N
overlap(:,ii) = abs(integrateMeshFunction(mesh.p,mesh.t,bsxfun(@times,eigstruct.rho(:,ii),eigstruct.phi(:,1:N))));
end
%display overlap (uncomment below)
%format short; disp(overlap)