function [sigma,sigmaxx,sigmaxy,Nmax] = conducMagnetoFull(ene,ef,eneB,gam,kbT,Nmax)
%Call:  [sigmaxx,sigmaxy,Nmax] = conducMagnetoFull(ene,ef,eneB,gam,kbT,Nmax)
%Calculates the local graphene conductivity in the presence of a finite
%magnetic B-field oriented perpendicular to the sheet.
%The unit of ene (hbar*omega), ef (Fermi level), gam (hbar*gamma), eneB
%(a characteristic energy/frequency of the B-field in graphene, given by
%eneB = hbar*omegaB = sqrt(2*e*B*hbar)*vf], and kBT (kb*T) should all be
%identical. Otherwise the unit is non-essential. The conductivity is output 
%in SI units. 

%Default value for Nmax based on maximum requested energy and eneB
if ~exist('Nmax','var') || isempty(Nmax)
    Nmax = ceil( 300*(max(ene) + abs(ef))^2/eneB^2 ); %The xx-part converges rather slowly, hence the large cutoff value
                                                     %Ideally, this issue should be fixed by an asymptotic scheme where 
                                                     %the sum for n>Nmax is converted to an integral over energy. [in practice this could be done as in Gusyin Eq. (13)] 
end


nlist = 0:1:(Nmax+1);
fintra = sqrt(nlist(1:end-1)+1)-sqrt(nlist(1:end-1));
finter = sqrt(nlist(1:end-1)+1)+sqrt(nlist(1:end-1));
En = sqrt(nlist)*eneB; %All the positive Landau levels from 0 to Nmax+1
fdp = fermidirac(En,ef,kbT); %Occupation of all positive landau levels (and zero)
fdn = fermidirac(-En,ef,kbT); %Occupation of all negative landau levels (and zero)
eneigam = ene+1i*gam;
eneigam2 = eneigam.^2;
eneB2 = eneB^2;
eneigam2_div_eneB2 = eneigam2/eneB2;

pref = 3.874045877403706e-05; %Prefactor of e^2/(2*pi*hbar) in SI units

sigmaxx = zeros(size(ene)); sigmaxy = zeros(size(ene)); %Preallocate
for n = 1:Nmax
    sigmaxx = sigmaxx + ...
              ( ( fdp(n) - fdp(n+1) ) + ( fdn(n+1) - fdn(n) ) ) ./ ...
              ( ( eneigam2 - fintra(n)^2*eneB2 ) * fintra(n) )  +  ...
              ( ( fdn(n) - fdp(n+1) ) + ( fdn(n+1) - fdp(n) ) ) ./ ...
              ( ( eneigam2 - finter(n)^2*eneB2 ) * finter(n) ); 
end
for n = 1:ceil(Nmax/50) %This part converges MUCH faster than the xx component 
	sigmaxy = sigmaxy + ...
              ( ( fdp(n) - fdp(n+1) ) - ( fdn(n+1) - fdn(n) ) ) .* ...
                (   1./(eneigam2_div_eneB2 - fintra(n)^2) + ...
                    1./(eneigam2_div_eneB2 - finter(n)^2)   );
end
sigmaxx = pref*1i*eneB*eneigam.*sigmaxx;
sigmaxy = pref*sigmaxy;
%The conductivity tensor
sigma(1,1,:) = sigmaxx; 
sigma(1,2,:) = sigmaxy; 
sigma(2,1,:) = -sigmaxy; %The sign is required by the Onsager relation sigma_xy(B) = sigma_yx(-B) = -sigma_xy(B), see e.g. Grosso & Paravicini or Ashcroft and Mermin (Eq. 13.90)
sigma(2,2,:) = sigmaxx; 

function fd = fermidirac(ene,ef,kbT)
%Fermi-Dirac distribution function at energy ene for a Fermi energy
%(chemical potential, really) ef, and at temperature kbT, all supplied in
%equal units.

if ~isempty(kbT) && kbT ~= 0
    fd = 1./(1+exp((ene-ef)/kbT));
else
    fd = double(ene<=ef);
end