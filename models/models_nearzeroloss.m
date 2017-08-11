function models = models_nearzeroloss(B)

ef_eV = 0.2; L = 400e-9;
models(1) = struct('ef_eV',ef_eV,'L',L,'gam_eV',1.5e-4 ,'ene_eV',60.3e-3,'B',B);

