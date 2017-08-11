function f = models_adiabaticAbsorb_bump_v2(p,ene_eV,ef_eV,B,gam_eV,gambump_eV)

% Adiabatic absorption bump configuration in xy-space
x0 = [-12.75,-1]; 
y0 = [-10.25,-3.5];
xcut = [1.5,3]; ycut = [1,1]; 
HeaviS = @(vals,val0,cut) (1 + erf(-(vals-val0)./cut))/2; 
Hx = HeaviS(p(:,1),x0(2),xcut(2)).*HeaviS(-p(:,1),-x0(1),xcut(1));  %xcut*(1+20*HeaviS(Y,y0(1),ycut(1)))
Hy = HeaviS(p(:,2),y0(2),ycut(2)).*HeaviS(-p(:,2),-y0(1),ycut(1)); 
H = Hx.*Hy; 

% Material & bump properties
gam_eV = gam_eV + gambump_eV*H; 
sigma0 = conducIntra(ene_eV,ef_eV);
sigma = conducMagnetoClassical(ene_eV,ef_eV,B,gam_eV);

%"Occupation" matrix at every vertex point
f = sigma./sigma0;