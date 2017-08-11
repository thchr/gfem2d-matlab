function [oscstrength_norm,zeta] = optim_oscstrength(periphery,polar,hdata,plotmesh)

%Default parameters
if ~exist('polar','var') || isempty(polar)
    polar = 1; %1 = x, 2 = y
end
if ~exist('hdata','var') || isempty(hdata)
    hdata = 5/size(periphery,1);
end
if ~exist('plotmesh','var') || isempty(plotmesh)
    plotmesh = 0;
end

%% MESH
%Calculate oscillator strengths
if isstruct(periphery)
    mesh = periphery; clear periphery; 
else
    [~,mesh.p,mesh.t] = evalc('geomPolygon(periphery,hdata)');
end

if plotmesh == 1
    %Maximize window via subaxis
    subaxis(1,1,1,'MT',.04,'MB',.07,'ML',.07,'MR',.04)
    %Draw mesh
    trimesh(mesh.t,mesh.p(:,1),mesh.p(:,2),'Color',[0.0784,0.0784,0.0863],'LineWidth',.2); %Plot with an off-black color [from flatcolors()]
    %Clean up
    axis equal; box on;
    xlim(minmax(mesh.p(:,1))+[-1,1]*max(abs(mesh.p(:)))/15)
    ylim(minmax(mesh.p(:,2))+[-1,1]*max(abs(mesh.p(:)))/15)
    set(gca,'Fontsize',8,'LineWidth',.2); drawnow;
end

%% SOLVE PROBLEM 

[~,eigstruct] = calcPolarizability(mesh); 

A = sum(meshArea(mesh.p,mesh.t));
oscstrength_norm = eigstruct.oscstrength(:,polar)/A;
disp(max(oscstrength_norm))
if nargout >= 2
    zeta = eigstruct.zeta;
end


