function [p,t,nodes,cnct] = geomFancyInclusion(bnd,incl,hdatafun,dispmesh)
%CALL: [p,t,nodes,cnct] = geomFancyInclusion(bnd,incl,hdatafun,dispmesh)
%DESCRIPTION: Meshes a region indicated by an exterior boundary 'bnd', and
%a number of inclusions - which may intersect with 'bnd' - that introduce
%holes or "bites in" in the geometry specified by 'bnd'.
%INPUT: bnd  | 2xN array which indicates the exterior boundary
%       incl | M element cell array of inclusions (each of which are 2xK
%              arrays)
%       hdatafun | (x,y) function specifying face size function ('const' or
%                  empty/unspecified for constant, default size function)
%       dispmesh | 0 (no) or 1 (yes), indicating whether to plot mesh
%OUTPUT : p & t | triangulated mesh variables
%         nodes | boundary variables used to create mesh
%         cnct  | connectivity array used with 'nodes'
%NOTE: 'bnd' and 'incl' (cell array of 2xN vectors) MUST be clockwise
%sorted andprovided in implicitly connected fashion (i.e., each neighboring
%element must connect cyclically)
%BY:                                     Thomas Christensen (23 May, 2016)


if ~exist('incl','var') || isempty(incl)
    bnd = [-1,-1; -1,1; 1,1; 1,-1];
end
if ~exist('incl','var') || isempty(incl)
    theta = linspace(0,2*pi,38); theta=theta(1:end-1).';
    incl{1} = [cos(theta),sin(theta)]/4;
    incl{end+1} = bsxfun(@plus,incl{1},[1,0]);
    incl{end+1} = bsxfun(@plus,incl{1},[-1,-1]);
    incl{end+1} = bsxfun(@plus,incl{1},[1,1]);
    incl{end+1} = bsxfun(@plus,incl{1},[0,1]);
end
if ~exist('dispmesh','var') || isempty(dispmesh)
    dispmesh = 1;
end

Nab = 10000; %Number of interpolation points between boundary nodes

%Find the intersections between the boundary and the inclusions
bndnew = []; flag = [];
for bb = 1:size(bnd,1)
    a = bnd(bb,:); b = circshift(bnd,-1,1);  b = b(bb,:);
    %partition the distance between boundary points into much finer regions,
    %allowing higher resolution detection of intersections
    a2b = [linspace(a(:,1),b(:,1),Nab).' , linspace(a(:,2),b(:,2),Nab).'];
    ain = 0; bin = 0; addb = []; adda = []; entry = []; exit = [];
    for ii = 1:numel(incl) %loop over inclusions
        inout = inpolygon(a2b(:,1),a2b(:,2),incl{ii}(:,1),incl{ii}(:,2)).';
        if any(diff(inout)==1)
            entry = [entry; find(diff(inout)==1)];
        end
        if any(diff(inout)==-1)
            exit = [exit; find(diff(inout)==-1)+1];
        end
        ain = max(ain,inout(1));
        bin = max(bin,inout(end));
    end
    if ain~=1 && (isempty(bndnew) || all(~ismember(bndnew,a,'rows')))
        entry = [1; entry];
    end
    if bin~=1 && (isempty(bndnew) || all(~ismember(bndnew,b,'rows')))
        exit = [exit; Nab];
    end
    %New boundary points
    bndadd = a2b(sort([entry;exit],'ascend'),:);
    bndnew = [bndnew;bndadd];
    
    %Indicate if one of the original boundary points lie inside an inclusion
    flag = [flag; zeros(size(bndadd,1),1)];
    if bin == 1; flag = [flag(1:end-1); mod(bb,size(bnd,1))+1];  end
end

%add midpoints to boundary, which will allow us to check whether a boundary
%point should connect to its neighbour or not
for bb = 1:size(bndnew,1)
    if flag(bb) == 0
        a = bndnew(bb,:); b = circshift(bndnew,-1,1);  b = b(bb,:);
        mid(bb,:) = (a + b)/2;
    else
        mid(bb,:) = bnd(flag(bb),:);
    end
end


%delete parts of inclusions outside original boundary
inclnew=[]; inclcnct = [];
for ii = 1:numel(incl)
    cncttemp = bsxfun(@plus,repmat( (1:size(incl{ii},1))',1,2),[0,1]); cncttemp(end,2) = 1;
    inout=inpolygon(incl{ii}(:,1),incl{ii}(:,2),bnd(:,1),bnd(:,2));
    kill = sort(find(~inout),'descend');
    for kk = 1:numel(kill); %Kill points outside bnd, and update counting scheme
        cncttemp(any(cncttemp==kill(kk),2),:) = [];
        cncttemp(cncttemp>kill(kk)) =  cncttemp(cncttemp>kill(kk))  - 1;
    end
    inclcnct = [inclcnct; cncttemp+size(inclnew,1)];
    inclnew  = [inclnew; incl{ii}(inout,:)];
end


%find boundary connections between newbnd and (proper parts of) inclusions
cnct=[];
for bb = 1:size(bndnew,1)
    for ii = 1:numel(incl)
        in = inpolygon(mid(bb,1),mid(bb,2),incl{ii}(:,1),incl{ii}(:,2)).';
        if in == 1; break; end
    end
    if bb < size(bndnew,1); bb2 = bb+1; else bb2 = 1; end
    if in == 0
        cnct = [cnct; bb,bb2];
    else
        [~,newcnct1] = min(pdist2(bndnew(bb,:),inclnew));
        [~,newcnct2] = min(pdist2(bndnew(bb2,:),inclnew));
        cnct = [cnct; bb,newcnct1+size(bndnew,1); newcnct2+size(bndnew,1),bb2];
    end
end

%add nodes and connections together
nodes = [bndnew; inclnew];
cnct = [cnct; size(bndnew,1)+inclcnct];

%There may be points which are too close in their intersection; set them equal,
%and then let checkgeometry() sort them out - the default tolerance is set
%by a fraction of the minimum distance between inclusions
minincldist = inf;
for cci = 1:size(inclcnct,1) %Finding the appropriate tolerance
    minincldist = min(minincldist,sqrt(diff(inclnew(inclcnct(cci,:),1))^2+...
        diff(inclnew(inclcnct(cci,:),2))^2));
end
for cc = 1:size(cnct,1);  %Finding points where there might be an issue
    cnctlen = sqrt(diff(nodes(cnct(cc,:),1))^2+diff(nodes(cnct(cc,:),2))^2);
    %If there is an issue we resolve by taking average value and equaling
    if cnctlen < minincldist/20;  %choice of tolerance
        newnode = sum(nodes(cnct(cc,:),:),1)/2;
        nodes(cnct(cc,:),:) = repmat(newnode,2,1);
    end
end


[nodes,cnct] = checkgeometry(nodes,cnct,[],[]); %now that unique() can work on the issue-elements, checkgeometry() can fix them


%plot mesh outline (depending on input)
if dispmesh == 1
    set_figsize([],50,25); hold on
    for cc = 1:size(cnct,1)
        plot(nodes(cnct(cc,:),1),nodes(cnct(cc,:),2),'-','color',[0.1173,0.2125,0.5275],'LineWidth',2.5)
    end;
    axis equal; box on;
    xlim(minmax(nodes(:,1))+[-1,1]*max(abs(nodes(:)))/15)
    ylim(minmax(nodes(:,2))+[-1,1]*max(abs(nodes(:)))/15)
    set(gca,'Fontsize',8,'LineWidth',.2); drawnow;
end

%default hdata function or specified
if ~exist('hdatafun','var') || isempty(hdatafun) || strcmpi(hdatafun,'const') %default (constant)
    hdata = struct('fun',@(x,y) ones(size(x))*minincldist*1.25);
else %specified
    hdata.fun = hdatafun; clear hdatafun;
end

%mesh the structure using mesh2d() and smooth with smoothmesh()
[p,t,stats] = mesh2d(nodes,cnct,hdata,struct('plot',false)); disp(stats)
[p,t] = smoothmesh(p,t);

%plot mesh (depending on input)
if dispmesh == 1;
    trimesh(t,p(:,1),p(:,2),'Color',[0.0784,0.0784,0.0863],'LineWidth',.2); % plot mesh
    plot(nodes(:,1),nodes(:,2),'o','Color',[0.1173,0.2125,0.5275],'MarkerSize',3);
    drawnow;
    hold off
end

%% old code used for testing and visualizing validity of implementation
% set_figsize([],20,20); hold on
% for cc = 1:size(cnct,1)
% 	plot(nodes(cnct(cc,:),1),nodes(cnct(cc,:),2),'-')
% end
%hold off