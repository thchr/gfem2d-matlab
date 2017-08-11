clear all; close all; clc;
%addpath(genpath('..'))
ConstantsUnits0;
cols = flatcolors; darkcol = cols{4}; cols(4) = [];

%% MESH
adivb = 29/14;
mesh = geomEllipse([1,adivb],[],.05,1);
[alpha,eigstruct] = calcPolarizability(mesh,20);

%% GRAPHENE
ef = .4;
kbt =kb_eV*300;
gam = 1*10^12/ef*hbar_eV;
ene = linspace(0,.8,15000);
sigmatot = conducLRA(ene,ef,gam,kbt);
normtermtot = real(sigmatot.^(-1)).^(-1);
absboundtot = normtermtot/c/eps0;
%%
Ll = [14,28]*1e-9;
colc = {cols{3},cols{13}};
for ll = 1:numel(Ll)
    L = Ll(ll);
    
    zeta = 2i*eps0*ene/hbar_eV*L./sigmatot;
    Q_abs = (ene/hbar_eV/c).*imag(alpha.x(zeta,L))/(pi*L^2*adivb);
    
    if ll == 1; set_figsize(2,24*.65,17*.65);
        plot(ene,absboundtot,'-','linewidth',1,'color',darkcol)
        hold on
        ylim([0,8])
        xlim(minmax(ene))
    end
    plot(ene,Q_abs,'Color',colc{ll},'LineWidth',1.5)
    if ll == numel(Ll)
        xlabel('Energy (eV)','Fontsize',10)
        ylabel('(\it\sigma_{\rmabs}/\itA\rm)_{max}','Fontsize',10)
        set(gca,'Fontsize',7,'LineWidth',.65,'Layer','Top')
    end
    
end
%export_fig('C:\Dropbox (MIT)\Projects\Limits in 2D plasmonics\AbsorptionBounds_14x29Ellipse','-pdf')