function [eneeig_eV,rhoeig,phieig,misc] = calcEigenEnergiesAny(mesh,models,bloch)
%CALL: [eneeig_eV,rhoeig,phieig,misc] = calcEigenEnergies(mesh,model,bloch)
%DESCRIPTION: Calculates eigenenergies (in eV) and eigendensities for a
%mesh (with field p and t) and a models (an array of structures) w/ fields:
%       models.type | 'isotropic' or 'magneto' or 'berry'
%       models.approx | 'intra' or 'full' (only if type = 'isotropic')
%       models.system | 'graphene' or '2deg' (only if type = 'magneto')
%       models.ef_eV | Fermi energy (in eV)
%       models.L | Characteristic length with which mesh is scaled (in m)
%       models.B | Applied magnetic field (in T) (only if type = 'magneto')
%       models.gam_eV | If type = 'isotropic' and 'approx' = 'intra', it is
%                       possible to supply an intrinsic loss rate in eV
%       models.gamc   | If models.type is 'magneto', it is possible to
%                       supply an intrinsic loss rate in units of the
%                       cyclotron frequency, i.e. gamc = gam_eV/omegac_eV
%       model.keep_eV | 2-element array, which indicates which eigenvalues/
%                       vectors to retain for output; specified as range in
%                       eV (i.e. as [min_eV,max_eV])
%       models.F | the dimensionless net Berry curvature (if type = 'berry)
%A 'bloch' structure must be passed if the calculation type exhibits any
%periodicity; if it does not, it can be empty. Its possible fields are:
%       bloch.type | '0D' for a finite calculation (default)
%                    '1D' for a ribbon calculation
%                    '2D' for a lattice calculation
%       bloch.k | k-vector [1x2 array] specifying the bloch momentum
%                 (necessary for '1D' and '2D')
%       bloch.n_genexpint | If bloch.type is '1D' this value must be passed
%                           to specify number of generalized exponential
%                           integral terms in the Taylor expansion of the
%                           long-range Coulomb term in the Ewald scheme
%       bloch.supercell | If bloch.type is '1D', we can alternatively
%                         calculate the Coulomb interaction by a supercell
%                         scheme; in that case, parse this value as a non-
%                         empty array with 2 elements indicating
%                         multiplication onto R{1} and R{2} (extension).
%
%WRITTEN: by Thomas Christensen                             04 May, 2016.

%Default values
if ~exist('bloch','var') || isempty(bloch)
    bloch.type = '0D';
end
if strcmpi(bloch.type,'1D') || strcmpi(bloch.type,'2D');  %If 1D or 2D, we explicitly include zero-energy solutions
    kzero = all(bloch.k == [0,0]);                        %if k == [0,0] (only relevant for isotropic calculation); there kzero = 1
else                                                      %If 0D, we don't encounter this scenario, and set kzero = 0
    kzero = 0;
end

%Write to prompt the purpose and setup considered
fprintf('Calculating eigenenergies of specified mesh (%g nodes and %g elements)\n',size(mesh.p,1),size(mesh.t,1));

%Check if any of the models require us to calculate DRa; then chckani = 1 (otherwise, if all are 'isotropic', chckani = 1)
for tt = {models(:).type}; if ~strcmpi(tt,'isotropic'); chckani = 1; break; else chckani = 0; end; end
%Construct the necessary matrices for the calculation
if chckani == 0 %Only isotropic models requested
    [V,DL,DRi] = calcSystemMatrices(mesh,bloch); DRa = [];
elseif chckani == 1 %At least one isotropic model requested
    [V,DL,DRi,DRa] = calcSystemMatrices(mesh,bloch);
end

%Preallocating;
eneeig_eV = cell(numel(models));
if nargout >= 2; rhoeig = cell(size(models)); end
if nargout >= 3; phieig = cell(size(models)); end
if nargout == 4; misc   = cell(size(models)); end

fprintf('Finding eigenvalues for each of the %g model-input(s)\n',numel(models))
%Looping over models
for mm = 1:numel(models)
    fprintf('\n   Model %g/%g | Parameters:\n',mm,numel(models));
    disp(models(mm)); %Print model structure as display [alternatively, try disp(struct2table(models(mm)))]
    %Switch based on model parameters
    if nargout == 1
        [eneeig_eV{mm}] = wrapModelSpecificCalc(V, DL, DRi, DRa, models(mm),kzero);
    elseif nargout == 2
        [eneeig_eV{mm},rhoeig{mm}] = wrapModelSpecificCalc(V, DL, DRi, DRa, models(mm),kzero);
    elseif nargout == 3
        [eneeig_eV{mm},rhoeig{mm},phieig{mm}] = wrapModelSpecificCalc(V, DL, DRi, DRa, models(mm),kzero);
    else
        [eneeig_eV{mm},rhoeig{mm},phieig{mm},misc{mm}] = wrapModelSpecificCalc(V, DL, DRi, DRa, models(mm),kzero);
    end
end

%If models is only a one-element array (i.e. only one structure) we don't
%output the results as cells
if numel(models) == 1
    eneeig_eV = eneeig_eV{1};
    if nargout >= 2;   rhoeig = rhoeig{1};  end
    if nargout >= 3;   phieig = phieig{1};  end
    if nargout >= 4;   misc = misc{1};      end
end

%%
%-------------------------------------------------------------------------
%----------------------------- SUB-FUNCTIONS -----------------------------
%-------------------------------------------------------------------------

%% ------------------------------------------------------------------------
function [eneeig_eV,rhoeig,phieig,misc] = wrapModelSpecificCalc(V, DL, DRi, DRa, model,kzero)
ConstantsUnits0;
switch model.type
    case 'isotropic' %Local isotropic calculations
        [zeta,rhoeig] = calcEigen(V,DL,DRi,[],[],kzero); %Retain zero solution if k = [0,0]
        
        switch model.approx
            case 'intra' %The intraband approximation
                if ~isfield(model,'gam_eV') || model.gam_eV == 0 %Calculation without intrinsic loss
                    eneeig_eV = sqrt( (e^2*model.ef_eV*ev2jo)/(2*pi*eps0)*(real(zeta)/model.L) ) / ev2jo;
                else %Calculation with intrinsic loss (via model.gam_eV)
                    eneeig_eV = sqrt( (e^2*model.ef_eV*ev2jo)/(2*pi*eps0)*(real(zeta)/model.L) / ev2jo^2 - (model.gam_eV/2)^2 ) ... %Real part
                        - 1i*model.gam_eV/2;                                                                                %Imag part
                end
            case 'full' %The full local approximation
                eneeig_eV = zeros(size(zeta)); %Preallocation
                if ~isfield(model,'kT_eV'); model.kT_eV = []; end
                guess = 0.1; %arbitrary starting guess
                for zz = 1:numel(zeta);
                    zetaf = @(ene_eV) 2*eps0*(ene_eV/hbar_eV)*model.L./imag(conducLRA(ene_eV,model.ef_eV,model.kT_eV));
                    eneeig_eV(zz) = fzero(@(ene_eV) real(zetaf(ene_eV) - real(zeta(zz))),guess);
                    guess = eneeig_eV(zz); %Update starting guess
                end
        end
        
        %Miscellaneous output
        if nargout == 4; misc.zeta = zeta; misc.model = model; end
        
    case 'magneto' %Intraband (local) anisotropic calculation for finite magnetic field (perpendicular)
        switch model.system
            case 'graphene'
                omegac = e*model.B*vf^2/(model.ef_eV*ev2jo); %Classical cyclotron: graphene
                sigmac = conducIntra(omegac*hbar_eV,model.ef_eV);    %Drude conductivity evaluated at cyclotron freq.
            case '2deg'
                omegac = e*model.B/model.meff;               %Classical cyclotron: 2DEG
                sigmac = 1i*eps0*(model.omegap_eV/hbar_eV)^2/omegac; %Drude conductivity evaluated at cyclotron freq.
        end
        zetac = 2i*eps0*omegac*model.L/sigmac;
        
        %Check if loss rate provided in eV (and not in gamc), if so, create
        %a gamc object
        if isfield(model,'gam_eV') && ~isfield(model,'gamc')
            model.gamc = model.gam_eV/(omegac*hbar_eV);
        end
        if nargout == 1
            lambda = wrapMagnetoEigen(V,DL,DRi,DRa,zetac,model);
        else
            [lambda,rhoeig] = wrapMagnetoEigen(V,DL,DRi,DRa,zetac,model);
        end
        eneeig_eV = lambda*omegac*hbar_eV;
        %Miscellaneous output
        if nargout == 4; 
            misc.zetac = zetac;
            misc.omegac_eV = omegac*hbar_eV;
            misc.lambda = lambda;
            misc.model = model;
        end
        
    case 'berry' %Local intraband treatment of finite net Berry flux (Sung/Rudner-style)
        sigmaf = conducIntra(model.ef_eV,model.ef_eV); 
        zetaf = 2i*eps0*(model.ef_eV/hbar_eV)*model.L/sigmaf;
        
        if nargout == 1
            lambda = wrapBerryEigen(V,DL,DRi,DRa,zetaf,model);
        else
            [lambda,rhoeig] = wrapBerryEigen(V,DL,DRi,DRa,zetaf,model);
        end
        eneeig_eV = lambda*model.ef_eV;
        %Miscellaneous output
        if nargout == 4; 
            misc.zetaf = zetaf; 
            misc.lambda = lambda;
            misc.model = model; 
        end
end

%Sort from small to large real part of the energy
[~,I]=sort(real(eneeig_eV)); eneeig_eV = eneeig_eV(I);
if isfield(model,'keep_eV')
    I2 = find(real(eneeig_eV)>=model.keep_eV(1) & real(eneeig_eV)<=model.keep_eV(2));
    eneeig_eV = eneeig_eV(I2);
end
%Save requested additional output
if nargout >= 2
    rhoeig = rhoeig(:,I);
    if isfield(model,'keep_eV'); rhoeig = rhoeig(:,I2); end;
end
if nargout >= 3; phieig = V*rhoeig; end

%% ------------------------------------------------------------------------
function [lambda,eigrho] = wrapMagnetoEigen(V,DL,DRi,DRa,zetac,model)

sz = size(V); %We need the size of V on several occasions; get it here

ticcompan = tic; %Bottom non-trivial row of the companion matrix
if ~isfield(model,'gamc') || model.gamc == 0 %Calculation without intrinsic loss
    brow = [(-1i/(2*pi*zetac))*(DL\(DRa*V))                   , (1/(2*pi*zetac))*(DL\(DRi*V)) + eye(sz)                 , zeros(sz)] ;
else                                         %Calculation with intrinsic loss (via model.gamc)
    brow = [(-1i/(2*pi*zetac))*(DL\(DRa*V - model.gamc*DRi*V)), (1/(2*pi*zetac))*(DL\(DRi*V)) + (1+model.gamc^2)*eye(sz), -2i*model.gamc*eye(sz) ];
end

%Companion matrix
C = [zeros(sz), eye(sz), zeros(sz); zeros(sz), zeros(sz), eye(sz); brow];
fprintf('   Companion matrix for cubic eigenvalue problem created in %.2f sec\n',toc(ticcompan))

%Solve the resulting 3N x 3N eigenvalue problem
fprintf('   Cubic eigenvalue problem of %g x %g system',sz(1),sz(2)); ticpeig = tic;
if nargout == 1 %Only compute eigenvectors if requested; takes a lot of additional time, which can otherwise be saved
    lambda = eig(C,'vector');
else
    [eigrho,lambda] = eig(C,'vector');
    eigrho = eigrho(1:sz(1),:); %Cut the lambda-repeated parts off (the last two thirds)
end

fprintf(' solved in ') %Print solve-time for eigenvalue problem
if toc(ticpeig) > 60; fprintf('%.1f min\n',toc(ticpeig)/60); else fprintf('%.1f sec\n',toc(ticpeig)); end;

%% ------------------------------------------------------------------------
function [lambda,eigrho] = wrapBerryEigen(V,DL,DRi,DRa,zetaf,model)

sz = size(V); %We need the size of V on several occasions; get it here

ticcompan = tic; %Bottom non-trivial row of the companion matrix
if ~isfield(model,'gamf') || model.gamf == 0 %Calculation without intrinsic loss
    brow = [(1/(2*pi*zetaf))*(DL\(DRi*V)), (-1i*model.F/(2*zetaf))*(DL\(DRa*V))] ;
else                                         %Calculation with intrinsic loss (via model.gamf = gamma_eV/ef_eV)
    error('Not yet implemented\n')
end

%Companion matrix
C = [zeros(sz), eye(sz); brow];
fprintf('   Companion matrix for quadratic eigenvalue problem created in %.2f sec\n',toc(ticcompan))

%Solve the resulting 2N x 2N eigenvalue problem
fprintf('   Quadratic eigenvalue problem of %g x %g system',sz(1),sz(2)); ticpeig = tic;
if nargout == 1 %Only compute eigenvectors if requested; takes a lot of additional time, which can otherwise be saved
    lambda = eig(C,'vector');
else
    [eigrho,lambda] = eig(C,'vector');
    eigrho = eigrho(1:sz(1),:); %Cut the lambda-repeated parts off (the last half)
end

fprintf(' solved in ') %Print solve-time for eigenvalue problem
if toc(ticpeig) > 60; fprintf('%.1f min\n',toc(ticpeig)/60); else fprintf('%.1f sec\n',toc(ticpeig)); end;