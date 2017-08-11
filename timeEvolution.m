function [T,rho,phiind] = timeEvolution(mesh,models,Tmax,phiext,dTphiext,init,bloch)
%Setup is limited to simple Drude description of graphene. Can be applied
%to general 2D Drude models of plasma freq omegap by using a Fermi energy
%ef = pi*eps0*hbar^2*omegap^2/e^2.

%% MISC "STARTUP"
ConstantsUnits0;

%Write to prompt the purpose and setup considered
fprintf(['Time-evolution via ode45 time-stepping requested:\n' ...
    '   Mesh with %g nodes and %g elements\n' ...
    '   %g model(s) considered\n' ...
    '   Time-span [%.1f,%.1f 1/omegaF] in %g recorded steps\n\n'],...
    size(mesh.p,1),size(mesh.t,1),numel(models),Tmax(1),Tmax(2),Tmax(3));

%Default values
if ~exist('bloch','var') || isempty(bloch)
    bloch.type = '0D';
end

%Central finite differences approximation of derivative of phiext at prespecified step size
if ~exist('dTphiext','var') || isempty(dTphiext)
    dT = 1e-6;
    dTphiext = @(T) (phiext(T+dT/2) - phiext(T-dT/2))/dT;
end

%Print-out function for ode45
%odeopts = ...; 

%% SYSTEM MATRICES
%Check if any of the models require us to calculate DRa; then chckani = 1 (otherwise, if all are 'isotropic', chckani = 1)
for tt = {models(:).type}; if ~strcmpi(tt,'isotropic'); chckani = 1; break; else chckani = 0; end; end
%Construct the necessary matrices for the calculation
if chckani == 0 %Only isotropic models requested
    [V,DL,DRi] = calcSystemMatrices(mesh,bloch);
elseif chckani == 1 %At least one isotropic model requested
    [V,DL,DRi,DRa] = calcSystemMatrices(mesh,bloch);
    Da = -full(DL\DRa);
end
Di = -full(DL\DRi); 
Nvert = size(Di,1);

%% LOOPING OVER MODELS
fprintf('\nPerforming numerical ODE time-stepping (via ode45) for %g model(s)\n',numel(models))
%Looping over models
for mm = 1:numel(models)
    model = models(mm);
    fprintf('\n   Model %g/%g | Parameters:\n',mm,numel(models));
    disp(model); %Print model structure as display [alternatively, try disp(struct2table(models(mm)))]
    
    gam_n = model.gam_eV/model.ef_eV;
    Q = e^2/(4*pi*eps0*hbar*model.L*(model.ef_eV/hbar_eV));
    %Construct the appropriate ode-function structure based on model parameters
    if strcmpi(model.type,'isotropic')
        DiV = Di*V;
        odefun = @(T,Y) odefunIsotropic(DiV,Di,gam_n,Q,phiext,T,Y);
        if ~exist('init','var') || isempty(init); init = zeros(Nvert*2,1); end
    elseif strcmpi(model.type,'magneto')
        omegac_n = e*model.B*vf^2/(model.ef_eV*ev2jo)*hbar_eV/model.ef_eV;
        [A,B,C,CV,D] = odecoefsMagneto(V,Di,Da,gam_n,Q,omegac_n);
        odefun = @(T,Y) odefunMagneto(A,B,C,CV,D,phiext,dTphiext,T,Y);
        if ~exist('init','var') || isempty(init); init = zeros(Nvert*3,1); end
    elseif strcmpi(model.type,'berry')
        error('Not implemented yet\n')
    else
        error('Model type is not supported (''models'' field ''type'' must be either ''isotropic'', ''magneto'', or ''berry'')\n')
    end
    
    [T{mm},rhoY] = ode45(odefun,linspace(Tmax(1),Tmax(2),Tmax(3)),init,odeset('OutputFcn',@odetpbar));
    rho{mm} = rhoY(:,1:Nvert).'; %Remove repetitions
    if nargout >= 3; phiind{mm} = V*rho{mm}/(4*pi*eps0); end
end

%don't need to output as cell array if there's only one model
if numel(models) == 1
    T = T{1}; rho = rho{1}; if nargout >= 3;  phiind = phiind{1}; end
end
%% FIRST-ORDER ODE MATRIX SETUP
function dY = odefunIsotropic(DiV,Di,gam_n,Q,phiext,T,Y)

Nvert = size(Di,1); dY = zeros(2*Nvert,1);

dY(1:Nvert) = Y(Nvert+1:end);
dY(Nvert+1:2*Nvert) = + (Q/pi)*DiV*Y(1:Nvert) - gam_n*Y(Nvert+1:end)  ...     %NOTENOTENOTE: It seems the sign in front of -Q*Di*V*Y has something to do with the divergent behavior!!!
                      + (Q/pi)*Di*phiext(T);


function [A,B,C,CV,D] = odecoefsMagneto(V,Di,Da,gam_n,Q,omegac_n)

A = 2*gam_n;
B = (omegac_n^2 + gam_n^2)*speye(size(Di)) - (Q/pi)*Di*V;
C = (Q/pi)*(-gam_n*Di + omegac_n*Da);
CV = C*V;
D = (Q/pi)*Di;

function dY = odefunMagneto(A,B,C,CV,D,phiext,dTphiext,T,Y)

Nvert = size(C,1); dY = zeros(3*Nvert,1);

dY(1:Nvert) = Y(Nvert+1:2*Nvert);
dY(Nvert+1:2*Nvert) = Y(2*Nvert+1:end);
dY(2*Nvert+1:end) = - CV*Y(1:Nvert) - B*Y(Nvert+1:2*Nvert) - A*Y(2*Nvert+1:end)  ...
    + D*dTphiext(T) - C*phiext(T);  %driving term
