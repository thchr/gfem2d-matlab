function models = models_tstep(B,gam_eV)

ef_eV = 0.2; L = 400e-9;
if B == 0
    models = struct('type','isotropic','ef_eV',ef_eV,'L',L,'gam_eV',gam_eV,'B',B);
else
    models = struct('type','magneto',  'ef_eV',ef_eV,'L',L,'gam_eV',gam_eV,'B',B);
end

