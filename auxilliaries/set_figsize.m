function hfig = set_figsize(fignum,width,height,opts)
%Initialize a figure window of number 'fignum' (set to [] if a new window
%without specified fignum is desired) with width and height specified by
%[width,height]. Additional options to fig may be set through opts. Figure
%color is set to white as standard.
%
%Syntax as follows:
%     hfig = = set_figsize(fignum,width,height,opts)
%width may be specified as string 'max' to indicate the figure should be
%maximized

if ~exist('width','var') || isempty(width) %Standard choice for figure width
    width = .9*8.64468;
end

if ~exist('height','var') || isempty(height) %Standard choice for figure width
    height = .7*8.64468;
end

if exist('fignum','var') && ~isempty(fignum)
    try close(fignum); end %Make sure the figure is closed
    hfig = figure(fignum);
else
    hfig = figure;
end

set(hfig,'color','white','units','centimeters')
if isnumeric(width) && isnumeric(height) %Specify figure dimensions
    set(hfig,'papersize',[width,height])
    pos = get(hfig,'position');
    set(hfig,'position',[27,17,width,height]);
    movegui(hfig,'northeast'); %Move figure to top-right corner
    drawnow; %Force a moving before any further command is executed
    
elseif strcmpi(width,'max') %Maximize figure window
    drawnow; %Command must be executed prior to maximizing figure (java complains otherwise)
    jFrame = get(handle(hfig), 'JavaFrame');
    jFrame.setMaximized(1);
end

if nargin == 4
    set(hfig,opts{:})
end
