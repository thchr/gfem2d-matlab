clear all; close all; clc; 
addpath(genpath('..'))
ConstantsUnits0; 
cols = flatcolors; darkcol = cols{4}; cols(4) = []; 
efL = 0.1:.1:0.5;
for ee = 1:numel(efL)
    ef = efL(ee); 
kbt =kb_eV*300;
gam = 1*10^12/ef*hbar_eV;
ene = linspace(0,1,100); 

Z0 = 1/c/eps0;

sigmatot = conducLRA(ene,ef,gam,kbt);
sigmaintra = conducIntra(ene,ef,gam,kbt);

normtermtot = real(sigmatot.^(-1)).^(-1); 
%normtermapprox = (real(sigmaintra.^(-1)) - real(sigmainter./sigmaintra.^2)).^(-1);
normtermlim = e^2*ef/(pi*hbar*gam); 
normtermapprox = (1/normtermlim + pi^2*hbar/4/e^2*(2*ene.^2*gam)/ef).^(-1);
normtermapprox = normtermlim + e^2/(4*pi*hbar*eps0*c)*gam/ef*(-3*(ene/gam).^2 +1)/Z0;

absboundtot = Z0*normtermtot; 
absboundapprox = Z0*normtermapprox;
absboundlim = Z0*normtermlim; 

if ee == 1; set_figsize(1,10,8.5); end
plot(ene,absboundtot,'-','linewidth',1.75,'color',cols{ee})
hold on 
plot(ene,absboundapprox,'-.','color',cols{ee}*.5+darkcol*.5,'linewidth',1.25)
plot(minmax(ene).*[1,.15],[1,1]*absboundlim,':','color',darkcol,'linewidth',.5)
ylim([0,12])
xlim(minmax(ene))
if ee == numel(efL)
xlabel('Energy (eV)','Fontsize',10)
ylabel('(\it\sigma_{\rmabs}/\itA\rm)_{max}','Fontsize',10)
set(gca,'Fontsize',7,'LineWidth',.65,'Layer','Top')
end
end
%export_fig('C:\Data\Dropbox (MIT)\Projects\Limits in 2D plasmonics\AbsorptionBounds_InterbandCorrection','-pdf')