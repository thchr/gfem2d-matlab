function f = models_adiabaticAbsorb_bump(p,ene_eV,ef_eV,B,gam_eV,gambump_eV)

% Adiabatic absorption bump configuration in xy-space
x0 = -8.75; 
y0 = [-3.25,2.75];
xcut = .65; ycut = [2,1.25]; 
HeaviS = @(vals,val0,cut) (1 + erf(-(vals-val0)./cut))/2; 
Hx = HeaviS(p(:,1),x0,xcut);  %xcut*(1+20*HeaviS(Y,y0(1),ycut(1)))
Hy = HeaviS(p(:,2),y0(2),ycut(2)).*HeaviS(-p(:,2),-y0(1),ycut(1)); 
H = Hx.*Hy; 

% Material & bump properties
gam_eV = gam_eV + gambump_eV*H; 
sigma0 = conducIntra(ene_eV,ef_eV);
sigma = conducMagnetoClassical(ene_eV,ef_eV,B,gam_eV);

%"Occupation" matrix at every vertex point
f = sigma./sigma0;