function sigmaRPA_T = conducRPA_Maldague(ene,q,ef,gam,kbT)
%CALL:        [sigma, chi_gam, chi_rpa] = conducRPA_Maldague(ene,q,gam,Ef)
%DESCRIPTION: Calculates the nonlocal, finite-temperature RPA conductivity 
%sigma(omega,q) of graphene for finite doping, decay rate, and nonzero 
%temperature.
%INPUT: 
%   ene | Optical frequency in eV (ene = hbar_eV*omega) [possibly an array]
%   q   | Crystal momentum in SI units (1/m)
%   ef  | Fermi energy in eV
%   gam | Decay rate in eV
%   kbT | Thermal energy kB*T in eV
%OUTPUT: 
%   sigma   | Conductivity in SI units
%METHOD: Uses an identity of Maldague's [see Ben Van Duppen et al., 
%2D Mater. 3, 015011, 2016 (doi.org/10.1088/2053-1583/3/1/015011), Eqs. 
%(26-29)] to get the temperature-dependent nonlocal conductivity of
%graphene from the zero-temperature nonlocal conductivity by a clever
%integral trick.
%LIMITATIONS: 
% -- So far, only ene can be array: all other input must be scalars!
% -- The function does not perform well for small wave vectors (q<0.001kf),
%    and may also have issues at very small temperatures (keep kbT>0.005ef
%    to be safe)
%IMPLEMENTATION: June, 2017                             Thomas Christensen

sigmaRPA_T = zeros(size(ene));

for ee = 1:numel(ene)
        fun = @(y) conducRPA(ene(ee),q,y,gam) ./ ... 
                  (4*kbT*cosh( (y - ef)/(2*kbT) ).^2 );

        sigmaRPA_T(ee) = integral(fun,-inf,inf); %Replace 'integral' by 'quadgk' if Matlab 
                                                 %version does not have new 'integral' routine
end
