function sigma = conducIntra(ene,ef,gam,kbT)
%CALL:    sigma = conducIntra(ene,ef,gam,kbT)
%DESCRIPTION: Calculates the intraband part of the local graphene 
%conductivity. The unit of ene (hbar*omega), ef (Fermi level), gam 
%(hbar*gamma), and kBT (kb*T) should all be  identical. Otherwise the unit 
%is non-essential. The conductivity is output in SI units. 

pref = 7.748091706447625e-05i; %Prefactor equals i*e^2/(pi*hbar) [SI]
ef = abs(ef); %pos. ef is assumed below, but results are the same for neg. ef

if nargin == 2
    sigma = pref*ef./ene;
elseif nargin == 3 
    sigma = pref*ef./(ene+1i*gam);
elseif nargin == 4
    if kbT ~= 0 || ef/(2*kbT) < 17
        sigma = 2*pref*kbT./(ene+1i*gam) .* log(2*cosh(ef/(2*kbT)));
    else %In this case, log(2*cosh(x)) = log(exp(x)) = x to within numerical accuracy
        sigma = pref*ef./(ene+1i*gam);
    end
end