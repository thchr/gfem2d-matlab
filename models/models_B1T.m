%Creates 20 models, which go through the 0 T, 0.5 T, 1 T, .., 10 T
%magneto-optic scenarios for the following setup parameters:

ef_eV = .2; L = 400e-9; 
omegac_eV = @(B) 5.403927588936477e-04*B/ef_eV;

clear models
models(1) = struct('type','magneto'  ,'ef_eV',ef_eV,'L',L,'B',1,'system','graphene','keep_eV',[0,ef_eV]);
