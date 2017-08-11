clc; clear all; close all; 

%% MISC 
%addpath(genpath('..')) % don't need this if all files are in same folder
addpath('auxilliaries') % need this if the above is not used (to include in path the folder 'auxilliaries')
cols = flatcolors; % some color choices
ConstantsUnits0; % various necessary constants

%% PARAMETER CHOICES
ef_eV = 0.2; % Fermi level (eV)
gam_eV = 1e-3; % intrinsic loss (eV)
kbT_eV = 0e-3; % temperature (eV)
B = 8; % magnetic field (Tesla)

ene_eV = linspace(0,1,150)*ef_eV; % driving frequency/energy (eV)

%% CONDUCTIVITY CALCULATION 

% calculate the conductivity in quantum and semiclassical response models
eneB_eV = sqrt(2*e*B*hbar)*vf/ev2jo; %The level splitting between 0th and 1st Landau level

[~, semiclas.sigmaxx, semiclas.sigmaxy, semiclas.enec_eV] = conducMagnetoClassical(ene_eV,ef_eV,B,gam_eV); % semiclassical 
[~, landau.sigmaxx  , landau.sigmaxy,   landau.summedLevels] = conducMagnetoFull(ene_eV,ef_eV,eneB_eV,gam_eV,kbT_eV); %quantum

%% PLOTTING
landau.lnsty = {'-','color',cols{14}};
semiclas.lnsty = {'--','color',cols{21}};
normfac = e^2/hbar; %Normalization (= sigma_0)

set_figsize(1,15,10)
for tt = 1:2 % loop over the 'xx' and 'xy' components of the conductivity
    if tt == 1;
        plotlandau = landau.sigmaxx;
        plotsemiclas = semiclas.sigmaxx;
        plottype = 'xx';
    elseif tt == 2
        plotlandau = landau.sigmaxy;
        plotsemiclas = semiclas.sigmaxy;
        plottype = 'xy';
    end
    % compute y-range for plotting
    yran.imag = minmax([minmax(imag(plotlandau)),minmax(imag(plotsemiclas))]);
    yran.imag = yran.imag + [-1,1]*max(abs(yran.imag))*.05;
    yran.real = minmax([minmax(real(plotlandau)),minmax(real(plotsemiclas))]);
    yran.real = yran.real + [-1,1]*max(abs(yran.real))*.05;
    
    % plot the imaginary part of the conductivity component 
    subaxis(2,2,1,tt,'SV',.025,'SH',.035,'MR',.015,'ML',.07);
    plot(ene_eV,imag(plotlandau)/normfac,landau.lnsty{:}); hold on
    plot(ene_eV,imag(plotsemiclas)/normfac,semiclas.lnsty{:}); hold off
    ylim(yran.imag/normfac)
    ylabel(['\sigma_{' plottype '}/\sigma_0'])
    if tt == 1;
        title('Imaginary part')
        set(gca,'XTickLabel',{''})
        legend({'Quantized response','Semiclassical response'},'Fontsize',8); legend boxoff
    else
        xlabel('Energy [eV]')
    end
    set(gca,'Fontsize',7,'LineWidth',.4,'TickLength',[.0125,.0125])
    
    % plot the real part of a conductivity component 
    subaxis(2,2,2,tt);
    plot(ene_eV,real(plotlandau)/normfac,landau.lnsty{:}); hold on
    plot(ene_eV,real(plotsemiclas)/normfac,semiclas.lnsty{:}); hold off
    ylim(yran.real/normfac)
    if tt == 1;
        title('Real part')
        set(gca,'XTickLabel',{''})
        text(ene_eV(end)-diff(ene_eV([1 end]))/25,(yran.real(end)-diff(yran.real([1 end]))/25)/normfac,...
            {['E_F = ' num2str(ef_eV) ' eV  |  T = ' num2str(kbT_eV*1e3) ' meV']; ...
             ['B = ' num2str(B) ' T  |  \gamma = ' num2str(gam_eV*1e3) ' meV']; ...
             },...
             'HorizontalAlignment','right','VerticalAlignment','Top','Fontsize',8)
    else
        xlabel('Energy [eV]')
    end
    set(gca,'Fontsize',7,'LineWidth',.4,'TickLength',[.0125,.0125])
end