function [sigma, chi_gam, chi_rpa] = conducRPA(ene,q,ef,gam)
%CALL:        [sigma, chi_gam, chi_rpa] = conducRPA(ene,q,gam,Ef)
%DESCRIPTION: Calculates the nonlocal, ZERO-temperature RPA conductivity 
%sigma(omega,q) of graphene for finite doping and decay rate.
%INPUT: 
%   ene | Optical frequency in eV (ene = hbar_eV*omega)
%   q   | Crystal momentum in SI units (1/m)
%   ef  | Fermi energy in eV
%   gam | Decay rate in eV
%OUTPUT: 
%   sigma   | Conductivity in SI units
%   chi_gam | Noninteracting density-density response function in SI units
%   chi_rpa | RPA-interacting density-density response function in SI units
%METHOD: The function implements the complex-frequency RPA result obtained 
%by Sernelius (the real-frequency results obtained by Stauber et al. and 
%Hwang & Das Sarma do not work in combination with the Mermin prescription) 
%and rigorously includes a nonzero loss rate via Mermin's relaxation-time
%approximation. [see B. E. Sernelius, Physical Review B 85, 195427 (2012)
%and N. D. Mermin, Phys. Rev. B 1, 2362 (1970)]
%
%IMPLEMENTATION: First implemented at some point in 2012; input/output
%updated in 2017                                    Thomas Christensen

%Constants
e = 1.60217657e-19;     %Elementary charge in Coulomb
ev2jo = 1.60217657e-19; %Conversion from eV to Joule by multiplication
hbar = 1.05457173e-34;  %Planck constant in Js
hbar_eV =  hbar/ev2jo;  %------||------- in eVs
aLC = .142e-9*sqrt(3);  %Lattice constant (1.42 Å);    vf depends on this value
tAB = 2.8;              %Hopping term fixed to 2.8 eV; vf depends on this value
vf = sqrt(3)*aLC*tAB/(2*hbar_eV); %Graphene Fermi velocity m/s

%Implementation assumes pos. ef, but the result is the same for neg. ef; 
%flip sign if necessary here.
ef = abs(ef); 

%Definitions
kf = ef/(hbar_eV*vf);
x = q./(2*kf); xp = 1./x; 
z = (ene+1i*gam)./(2*ef);

%Calculating pieces of chi_0 with nonzero tau
f = asin( (1-z).*xp ) + asin( (1+z).*xp ) ...
    - (z-1).*xp.*sqrt(1-( (z-1).*xp ).^2) ...
    + (z+1).*xp.*sqrt(1-( (z+1).*xp ).^2);

D0 = 4*ef*ev2jo/(2*pi*vf^2)/hbar.^2; %Density of states prefactor

%Noninteracting polarizability at complex frequency
chi_omegagam = -D0.*(1 + x.^2./(4*sqrt(x.^2-z.^2)).*(pi-f));

%Noninteracting polarizability at zero frequency
chi_omega0 = q./(2*pi*hbar*vf).*(2*xp - ...
    double((1-xp)>0).*(sqrt(1-xp.^2)-acos(xp)));

%Mermin correction for finite loss
chi_gam = (1+1i*gam./ene).*chi_omegagam ./ ...
    (1+1i*gam./ene.*chi_omegagam./chi_omega0);

%Associated noninteracting RPA conductivity
sigma = 1i*e^2*(ene/hbar_eV)./q.^2.*chi_gam;

%Associated interacting RPA polarizability, if requested (real part = 0
%indicates plasmon dispersion)
if nargout >= 3
    eps0 = 8.85418782e-12; %Vacuum permittivity in units of F/m
    Vq = e^2./(2*eps0*q);  %2D Coulomb interaction in momentum space
    chi_rpa = chi_gam./(1-Vq.*chi_gam);
end