function [sigma, sigmaxx, sigmaxy, enec_eV] = conducMagneto2DEG(ene_eV,omegap_eV,B,meff,gam_eV)
%Call:  [sigma, sigmaxx, sigmaxy, eneC_eV] = 
%                   calcMagnetoClassical(ene_eV,omegap_eV,B,meff,gam_eV)
%Computes the semi-classical intraband one-transition approximation of the
%classical dynamic Drude-Hall response of a 2D parabolic electron gas of
%electronic charge carriers with effective mass meff. Additional inputs are
%'ene_eV' (hbar*omega), plasma-frequency 'omegap_eV' (hbar*omegap), and 
%loss-rate 'gam_eV' (hbar/tau), all supplied in units of eV. The B-field is
%supplied in Tesla and the mass in kg (SI-units). Zero temperature (T=0) 
%is assumed. The anisotropic conductivity 
%       sigma = [sigmaxx,sigmaxy;-sigmaxy,sigmaxx] 
%is output, with its individual elements accesible also as second and third
%elements. Finally, the characteristic cycltron frequency (in eV) enec_eV = 
%hbar*omegac = e*B*hbar/meff is output.

ConstantsUnits0; %Units

omegac = e*B/meff; %Classical cyclotron frequency for mass = meff
enec_eV = hbar_eV*omegac; %Cyclotron energy
varomega = enec_eV./(ene_eV + 1i*gam_eV); %Cyclotron-normalized (inverse) frequencies

pref = 1i*eps0*(omegap_eV/hbar_eV)^2/omegac; %Prefactor (usual lossless Drude conductivity at omega=omegac) in SI units
common = pref*varomega./(1-varomega.^2);  %Common product in sigmaxx and sigmaxy

%Drude and Hall components of the frequency
sigmaxx = common;
sigmaxy = common.*(-1i*varomega);

%Conductivity tensor
sigma(1,1,:) = sigmaxx; 
sigma(1,2,:) = sigmaxy;  
sigma(2,1,:) = -sigmaxy; %Due to Onsager relation; see e.g. Ashcroft & Mermin
sigma(2,2,:) = sigmaxx;  