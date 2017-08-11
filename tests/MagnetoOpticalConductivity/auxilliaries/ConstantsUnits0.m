%-----DECLRARE AS GLOBAL VARIABLES-----
%global tAB_eV hbar hbar_eV ev2jo e c eps0 vf kb kb_eV


%-----CONSTANTS-----
aLC = .142e-9*sqrt(3); %Lattice constant i units of nm
tAB_eV = 2.8; %We fix hopping term to 2.8 eV and find vf from this value
ev2jo = 1.60217657e-19; %1 eV in Joule
hbar = 1.05457173e-34; %Planck constant in Js
hbar_eV =  hbar/ev2jo; %------||------- in eVs
kb = 1.3806488e-23; %Boltzmann constant in J/K
kb_eV = kb/ev2jo;   %--------||-------- in eV/K
e = 1.60217657e-19; %Elementary charge in Coulomb
me = 9.10938356e-31; %Electron rest mass in kg
c = 2.99792458e8; %Speed of light in m/s
lambda2jo = @(lambda) 2*pi*hbar*c./lambda; %Function that converts from lambda [in m] to energy [in J]
lambda2eV = @(lambda) 2*pi*hbar_eV*c./lambda; %Function that converts from lambda [in m] to energy [in J]
eps0 = 8.85418782e-12; %Vacuum permittivity in units of F/m
selfinteraction = 0.58*4.35974417e-18; %Units of Joule (corresponding to Coulomb self-interaction energy)
CoulombPref = e^2/(4*pi*eps0); %Units of C^2/[F/m] = Jm
A_atom = sqrt(3)/4*aLC^2;%Area occupied per carbon atom in hexagonal lattice
vf = sqrt(3)*aLC*tAB_eV/(2*hbar_eV); %Graphene Fermi velocity m/s
eV2THz = 1/hbar_eV/(2*pi)/1e12; %Converts, upon multiplication, from eV to THz
a0 = 4*pi*eps0*hbar^2/me/e^2; %Bohr radius in m
