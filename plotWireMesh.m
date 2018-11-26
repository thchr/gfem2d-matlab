function plotWireMesh(p,t,fignum,nonewfig)
%CALL:       plotWireMesh(p,t,fignum,nonewfig)
%DESCRIPTION: Displays a mesh of vertices (p) and connections (t) in a 
%nice-looking plot

if ~exist('nonewfig','var') || nonewfig ~= 1
    %Check input
    if ~exist('fignum','var'); fignum = []; end
    %Open figure window
    set_figsize(fignum,25,25);
    %Maximize window via subaxis
    subaxis(1,1,1,'MT',.04,'MB',.07,'ML',.07,'MR',.04)
end

%Draw mesh
trimesh(t,p(:,1),p(:,2),'Color',[0.0784,0.0784,0.0863],'LineWidth',.2); %Plot with an off-black color [from flatcolors()]

if ~exist('nonewfig','var') || nonewfig ~= 1 %Only if we're opening a new figure
%Clean up
axis equal; box on;
xlim(minmax(p(:,1))+[-1,1]*max(abs(p(:)))/15)
ylim(minmax(p(:,2))+[-1,1]*max(abs(p(:)))/15)

set(gca,'Fontsize',8,'LineWidth',.2)
end
drawnow;