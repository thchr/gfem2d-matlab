function models = models_adiabaticAbsorb_v3(B)

% The necessary, trivial parameters for the basic setup
ef_eV = 0.2; L = 400e-9; ene_eV = 60.3e-3; gam_eV = 1e-4; 
sigma0 = conducIntra(ene_eV,ef_eV,gam_eV);

basicmodel = struct('ef_eV',ef_eV,'L',L,'ene_eV',ene_eV,'gam_eV',gam_eV,...
                    'sigma0',sigma0,'B',B,'gambump_eV',[],'f',[]);

%The bump maximum loss-rate
gambump_eV = 1e-3; 

%Creating a function handle for the occupation function
for mm = 1:numel(gambump_eV); 
    models(mm) = basicmodel;
    models(mm).gambump_eV = gambump_eV(mm);
    models(mm).f = @(p) models_adiabaticAbsorb_bump_v3(p,...   
                            models(mm).ene_eV,...
                            models(mm).ef_eV,...
                            models(mm).B,...
                            models(mm).gam_eV,...
                            models(mm).gambump_eV);
end

