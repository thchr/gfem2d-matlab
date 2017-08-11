function models = models_loss(gam_meV)
%Creates 1D models, which go through 0,4, & 8 T in B, and also includes
%finite loss, as specified in meV (note the milli-eV choice!!)

ef_eV = .2; L = 400e-9; gam_eV = gam_meV*1e-3; 
B = [0,4,8]; 
omegac_eV = 5.403927588936477e-04*B/ef_eV;
gamc = gam_eV./omegac_eV; 
clear models

models(1) = struct('type','isotropic','ef_eV',ef_eV,'L',L,'B',B(1),'gam_eV',gam_eV,'gamc',[], ...
                   'system',[]        ,'approx','intra','keep_eV',[0,.6*ef_eV]); %The lower cutoff is the cyclotron frequency
models(2) = struct('type','magneto'  ,'ef_eV',ef_eV,'L',L,'B',B(2),'gam_eV',[]    ,'gamc',gamc(2),...
                   'system','graphene','approx',[]     ,'keep_eV',[0,.6*ef_eV]); 
models(3) = struct('type','magneto'  ,'ef_eV',ef_eV,'L',L,'B',B(3),'gam_eV',[]    ,'gamc',gamc(3),...
                   'system','graphene','approx',[]     ,'keep_eV',[0,.6*ef_eV]); %---||---