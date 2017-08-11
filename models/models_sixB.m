%Creates 20 models, which go through the 0 T, 0.5 T, 1 T, .., 10 T
%magneto-optic scenarios for the following setup parameters:

ef_eV = .2; L = 400e-9; 
omegac_eV = @(B) 5.403927588936477e-04*B/ef_eV;

clear models
models(1) =      struct('type','isotropic','ef_eV',ef_eV,'L',L,'B',0,'system',[]        ,'approx','intra','keep_eV',[0,ef_eV]);

magnetodefault = struct('type','magneto'  ,'ef_eV',ef_eV,'L',L,'B',4,'system','graphene','approx',[]     ,'keep_eV',[0,ef_eV]);
Blist = 2:2:10;
for bb = 1:numel(Blist);
    models(end+1) = magnetodefault; 
    models(end).B = Blist(bb); 
end

clear magnetodefault bb Blist