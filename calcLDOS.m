function [LDOS,LDOSnr,LDOSr] = calcLDOS(mesh,r,pol,models)
%CALL: [LDOS,LDOSnr,LDOSr] = calcLDOS(mesh,r,pol,ene,model)

ConstantsUnits0;

% DEFAULTS AND VARIABLE-CHECK
if ~exist('r','var') || isempty(r)
    r = [1,2,.2; 3,4,.3; 5,6,.4; 0,0,.5];
end
if ~exist('polar','var') || isempty(pol)
    pol = [1,0,0;0,1,0;0,0,1];
end
numpol = size(pol,1);
if ~exist('models','var') || isempty(models)
    models.ef = 0.2; models.gam = 1e-3; models.kbT = kb_eV*300; 
    models.ene = linspace(0.05,1,50)*models.ef; models.L = 400e-9; 
    models.sigma = conducLRA(models.ene,models.ef,models.gam,models.kbT); 
end

% EXTERNAL POTENTIAL 
phiext = zeros(size(mesh.p,1),numpol*size(r,1)); %preallocate
for rr = 1:size(r,1)
    for pp = 1:size(pol,1)
        phiext(:,(rr-1)*numpol+pp) = calcDipoleSource([mesh.p,zeros(size(mesh.p,1),1)],r(rr,:),pol(pp,:));
    end %Ordering is [(r1,pol1),(r1,pol2),..,(r1,polN) , (r2,pol1),(r2,pol2),..,(r2,polN), ... , (rM,pol1),(rM,pol2),..,(rM,polN)]
end     %With phiext size equaling P x N*M, for P mesh nodes, N polarizations, and M source points

% SYSTEM MATRICES
V = calcCoulomb(mesh.p,mesh.t);
Vmap = calcCoulombGeneralPoints(mesh.p,mesh.t,r); 
[DR,DL] = calcDifferential(mesh.p,mesh.t);

% RESPONSE
%Compute response due to each dipole in each model, and find induced field
for mm = 1:numel(models)
    for ee = 1:numel(models(mm).ene)
        if isfield(models(mm),'anisotropy') && models(mm).anisotropy == 1;
            [~,DR] = evalc('calcDifferential(mesh.p,mesh.t,models(mm).f(:,:,ee))'); %Run without output to command line
        end

        rho =  ( (DL+1i*models(mm).sigma(ee)/(4*pi*eps0*models(mm).L*models(mm).ene(ee)/hbar_eV)*DR*V) ) \ ... %We write it this (lengthy) way, to avoid potential memory issues of the mesh is very large.
               ( -1i*models(mm).sigma(ee)/(models(mm).ene(ee)/hbar_eV*models(mm).L^2)*DR*phiext ); 
        phiind = Vmap*rho;  %THIS NEEDS FIXING (to match r,pol-ordering)
        %ACTUALLY, WHAT YOU NEED IS THE ELECTRIC FIELD, NOT THE POTENTIAL
    end
end

