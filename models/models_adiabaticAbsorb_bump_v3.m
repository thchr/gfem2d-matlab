function f = models_adiabaticAbsorb_bump_v3(p,ene_eV,ef_eV,B,gam_eV0,gambump_eV)

% Adiabatic absorption bump configuration in xy-space
x0 = [-9.5,7]; 
y0 = [-10.25,-3.75];
xcut = [1.5,3]; ycut = [1,1]; 
HeaviS = @(vals,val0,cut) (1 + erf(-(vals-val0)./cut))/2; 
Hx = HeaviS(p(:,1),x0(2),xcut(2)).*HeaviS(-p(:,1),-x0(1),xcut(1));  %xcut*(1+20*HeaviS(Y,y0(1),ycut(1)))
Hy = HeaviS(p(:,2),y0(2),ycut(2)).*HeaviS(-p(:,2),-y0(1),ycut(1)); 
H = Hx.*Hy; 

% Material & bump properties
gam_eV = gam_eV0 + gambump_eV*H; 
sigma0 = conducIntra(ene_eV,ef_eV,gam_eV0);
sigma = conducMagnetoClassical(ene_eV,ef_eV,B,gam_eV);

%"Occupation" matrix at every vertex point
f = sigma./sigma0;