function sigma = conducInter(ene,ef,gam,kbT)
%CALL:    sigma = conducInter(ene,ef,gam,kbT)
%DESCRIPTION: Calculates the interband part of the local graphene
%conductivity. The unit of ene (hbar*omega), ef (Fermi level), gam 
%(hbar*gamma), and kBT (kb*T) should all be  identical. Otherwise the unit
%is non-essential. The conductivity is output in SI units. 


pref = 7.748091706447625e-05; %Prefactor equals e^2/(pi*hbar) [SI]
ef = abs(ef); %pos. ef is assumed below, but results are the same for neg. ef

if nargin == 2 || (gam == 0 && kbT == 0) %Lossless response at zero-temperature (gamma,kbT = 0)
    twoef = 2*ef;
    sigma = pref*0.25*( pi*double(ene>twoef) + 1i*log( abs((twoef-ene)./(twoef+ene)) ) );
    
elseif nargin == 3 || kbT == 0 %Lossy response at zero temperature (gamma ~= 0, kbT = 0)
    twoef = 2*ef; 
    sgn = double( -(ene>twoef) + (ene<=twoef) );  %sign of 2*ef-ene
    sigma = pi * double( ene>twoef ); %Interband Landau term (does NOT account for loss; this is consistent with "Hanson"-treatment in finite-T case below)          
    sigma = sigma - ... %Add the analytically derived finite-loss term
                 1i * sgn .* log( ( ene + 1i*gam + sgn.*twoef )./ ...
                                  ( -sgn.*(ene + 1i*gam) + twoef ) ); 
    sigma = pref/4*sigma; %Restore remaining prefactors
    
    %This implementation is so-chosen to be CONSISTENT with the
    %finite-temperature, finite-loss [(gam,kbT) ~= 0] case --- in fact, it 
    %is derived from that formula (see below). It has been _carefully_
    %checked that the two yield _identical_ results in the low-T limit.
    %Don't change any of the boolean operators; there is a subtlety at 
    %ene = 2*ef, and the current implementation ensures that no i*pi terms
    %give rise to issues (i.e. they cancel appropriately in this way).
    %Implemented on December 9, 2016.
    
elseif nargin == 4 %Account for both finite temperature and loss   
    integrand = @(enep,enepp) ( Hfunc(enepp,ef,kbT) - Hfunc(enep/2,ef,kbT) ) ./ ( (enep+1i*gam).^2 - 4*enepp.^2 );
    sigma_f = @(enep)  pref * ( pi/4*Hfunc(enep/2,ef,kbT) + 1i*(enep + 1i*gam).*quadgk(@(enepp) integrand(enep,enepp), 0, Inf) );
    
    sigma = zeros(size(ene));
    for ee = 1:length(ene)
        sigma(ee) = sigma_f(ene(ee));
    end
    
    %We are adding phenomenological loss to the inter-band (integral) term  %%% NOT the "bump"/step'ish-function term!
    %following the prescription by G.W. Hanson in 'Quasi-transverse 
    %electromagnetic modes supported by a graphene parallel-plate
    %waveguide', Jour. of Appl. Phys. 104, 083414 (2008).
end


function H = Hfunc(ene,ef,kbT)
%Population difference between energies -ene and +ene. Works also for T=0.
%Numerically stable by introducing a cutoff.

cuttol = 710;
EkT = ene./(kbT);
if kbT ~= 0
    
    H(EkT<cuttol) = sinh(EkT(EkT<cuttol))./(cosh(ef./(kbT))+cosh(EkT(EkT<cuttol)));
    
    H(EkT>cuttol) = 1-1./(exp((ene(EkT>cuttol)-ef)./(kbT))+1); %If EkT > cuttol then cosh, sinh
                                            %will fail and give Inf or NaN - however, in this 
                                            %limit the value of H is simply 1-n[F](E)
    
    H(isnan(H)) = 1; %If the returned answer is NaN this indicates that cosh
                     %could not be evaluated due to overflow - i.e., if only
                     %we assume that Ef is not absurdly large or kB*T absurdly
                     %small compared to E or Ef, then this means that H = 1.
                     %To emphasize: This approximation breaks down at very
                     %low temperatures, e.g. on the order of kB*T << 0.01 Ef.
                     
    if(sum(isnan(H)) ~= 0)
        disp('some values were NaN - corrected to be 1')
    end
else %In the limit T->0 then we have H(E) = 1 - n[F](E) with n[F](E) being 
     %the Fermi Dirac function, i.e. a step function; theta(E-E[F])
    H(ene<ef) = 0;
    H(ene>=ef) = 1;
end