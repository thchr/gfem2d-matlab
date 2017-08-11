function [sigma, sigmaxx, sigmaxy, enec_eV] = conducMagnetoClassical(ene_eV,ef_eV,B,gam_eV)
%Call:  [sigma, sigmaxx, sigmaxy, eneC_eV] = 
%                              calcMagnetoClassical(ene_eV,ef_eV,B,gam_eV)
%Computes the semi-classical intraband one-transition approximation of the
%full magneto-optical conductivity of graphene for energies 'ene_eV'
%(hbar*frequency), Fermi-energy 'ef_eV', and loss-rate 'gam_eV' (hbar/tau), 
%all supplied in units of eV. The B-field is supplied in Tesla (SI-unit).
%The temperature is assumed vanishing (T=0).
%The anisotropic conductivity sigma = [sigmaxx,sigmaxy;-sigmaxy,sigmaxx] is
%output, with its individual elements accesible also as second and third
%elements. Finally, the characteristic cycltron-like semi-classical
%magnetic energy enec_eV = hbar*omegac = e*B*vf^2*hbar/ef is output in eV.
%See e.g. Wang, Apell, & Kinaret PRB 86, 215450 (2012) Eqs. (15).


ConstantsUnits0; %Units

omegac = e*B*vf^2/(ef_eV*ev2jo); %Classical cyclotron frequency for mass = ef/vf^2
enec_eV = hbar_eV*omegac; %Cyclotron energy
varomega = enec_eV./(ene_eV + 1i*gam_eV); %Cyclotron-normalized (inverse) frequencies

pref = (1i*e^2*ef_eV*ev2jo)/(pi*hbar^2*omegac); %Prefactor (usual lossless intraband conductivity at omega=omegac) in SI units
common = pref*varomega./(1-varomega.^2); %Common product in sigmaxx and sigmaxy

%Drude and Hall components of the frequency
sigmaxx = common;
sigmaxy = common.*(-1i*varomega);

%Conductivity tensor
sigma(1,1,:) = sigmaxx; 
sigma(1,2,:) = sigmaxy;  
sigma(2,1,:) = -sigmaxy; %Due to Onsager relation; see e.g. Ashcroft & Mermin
sigma(2,2,:) = sigmaxx;  