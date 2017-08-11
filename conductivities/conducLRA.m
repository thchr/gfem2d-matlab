function [sigma,sigmaIntra,sigmaInter] = conducLRA(ene,ef,gam,kbT)
%CALL:      [sigma,sigmaIntra,sigmaInter] = conducLRA(ene,ef,gam,kbT)
%DESCRIPTION: Calculates the local graphene conductivity.
%The unit of ene (hbar*omega), ef (Fermi level), gam (hbar*gamma), and kBT
%(kb*T) should all be  identical. Otherwise the unit is non-essential. The
%conductivity is output in SI units. 

ef = abs(ef); %the conductivity results are the same for pos/neg ef - our 
              %implementation assumes positive ef, so we ensure correctness
              %here.
              
if nargin == 2 
    sigmaIntra = conducIntra(ene,ef); 
    sigmaInter = conducInter(ene,ef); 
elseif nargin == 3
    sigmaIntra = conducIntra(ene,ef,gam); 
    sigmaInter = conducInter(ene,ef,gam); 
elseif nargin == 4
    sigmaIntra = conducIntra(ene,ef,gam,kbT); 
    sigmaInter = conducInter(ene,ef,gam,kbT); 
else 
    error('At least two arguments are required for conducLRA')
end

sigma = sigmaIntra + sigmaInter;