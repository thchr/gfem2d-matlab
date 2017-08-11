function [blochmesh,origmesh,boundary] = geomRibbonInclusionScissor(R,rep2,inclusion,cut2,Nuc,hdatafun,showmesh)
%CALL: [blochmesh,origmesh] =
%      geomRibbonInclusionScissor(R,rep2,inclusion,cut2,Nuc,hdatafun,showmesh)
%DESCRIPTION: Creates a folded ring-like mesh for calculations on 1D
%             Bravais lattices, as suitable for a ribbon-like restriction
%             of a 2D Bravais lattice described by lattice vectors R. The
%             functions calcRibbonCoulomb and calcLatticeDifferential
%             require as input the blochmesh structure generated herein.
%NOTE DIFFERENCE: This function differs from geomRibbonInclusion() by
%             allowing the user to "cut off" parts of the lattice in the
%             R{2} direction; however, it is no longer possible to extend
%             along R{1}, and repetitions along R{2} must exceed unity. The
%             functionality is useful in engineering edge mode dispersion.
%INPUT: R | Two lattice vectors as cell array of 2D row vectors [default:
%           square lattice with lattice constant 1]
%       rep2 | integer which indicates the number of times each unit cell
%              is repeated in the ribbon-construction along R{2}. Minimum 
%              value is two. By default, we assume unit-repetition along 
%              R{1}, which indicates the direction of infinite extension.
%       inclusion | [N x 2] column vector of points inside unit cell which
%                   indicate an inclusion [default: 1/4 radius disk]. 
%                   [If inclusion is a disk, it is preferable that it has
%                   points at phi = 0 and = pi (use even # of points)]
%       cut2 | Scalar value in range [0,1] which indicates a percentage to
%              "cut" away from the top and bottom in units of R{2}.
%       Nuc | Total number of points along the unit cell [default: 250]
%       hdatafun | Meshing characteristics, supplied as function or string
%                  ('const' or 'antidot') [default = 'const']
%       showmesh | Set to 0 to avoid display of mesh, can otherwise be any
%                  of the strings 'blochmesh', 'remesh', 'origmesh', and
%                  'pairs'. [default = 'remesh']
%OUTPUT: blochmesh | The folded mesh supplied as a structure with fields
%                    p, t, area, remesh, and Rrib
%        origmesh | The original mesh (structure w. fields p, t, area).
%        boundary | Structure with elements 'boundary' and 'cnct' that
%                   allow for easy plotting of the ribbon unit cell.
%BY:                                   Thomas Christensen (19 May, 2016)

%% DEFAULT VALUES AND INITIALIZATION

%Default lattice vectors (if R is empty): a square of side 1
if ~exist('R','var') || isempty(R) || all(strcmpi(R,'square'))
    clear R;   R{1} = [1,0]; R{2} = [0,1];
elseif strcmpi(R,'triangular');
    clear R;   R{1} = [1,0]; R{2} = [cosd(60),sind(60)];
elseif all(size(R) == [1,1]) %If R is a scalar, we interpret it as the
    Rl = R; R=[];            %width/height of a square unit cell
    R{1} = Rl*[1,0]; R{2} = Rl*[0,1];
end
rotR{1} = [-R{1}(2),R{1}(1)]; %Rotated values for wiggling
rotR{2} = [-R{2}(2),R{2}(1)];

%Default R2 repetition value
if ~exist('rep2','var') || isempty(rep2)
    rep2 = 2;
elseif rep2 == 1
    error('Currently, the function does not handle rep2 = 1 as input, minimum is 2')
end
%Default cut2 value
if ~exist('cut2','var') || isempty(cut2)
    cut2 = 0;
end
rep = [1,rep2-cut2*2];

%Default total number of points along unit cell
if ~exist('Nuc','var') || isempty(Nuc)
    Nuc = 70;
end

%Number of points on each side from Nuc
ns = round([norm(R{1},2), norm(R{2},2)]/(2*(norm(R{1},2)+norm(R{2},2)))*Nuc);
rns = round(rep.*ns); %Repeated version, which accounts for the number of ribbon cells
for nn = 1:2   %Ensure that the number is odd
    if mod(rns(nn),2) == 0
        rns(nn) = rns(nn)+1;
    end
end


%Default inclusion and checking of duplicates in inclusion
if ~exist('inclusion','var') || strcmpi(inclusion,'disk') %Default inclusion
    ncirc = 36;                                           %is a 1/4 radius disk
    theta = linspace(0,2*pi,ncirc+1); theta = theta(1:end-1);
    inclusion = .25*[cos(theta);sin(theta)].';
else                                    %Ensure no duplicates in inclusion
    newinclusion = unique(inclusion,'rows','stable');
    if numel(newinclusion) ~= numel(inclusion);  %Warn if duplicates found
        inclusion = newinclusion;
        warning('Duplicates in ''inclusion'' were removed')
    else
        clear newinclusion
    end
end

fac = 1500; %Wiggling factor

%% DEFINE THE BOUNDARY OF THE RIBBON UNIT CELL (COLLISION DETECTION PRIMARILY)

%Boundary of the ribbon unit cell, which consist of several lattice unit cells (first x-coords, then y-coords)
bnd(:,1) = [linspace(0,rep(1)*R{1}(1),rns(1)).'; ...                     %"bottom"
    linspace(rep(1)*R{1}(1),rep(1)*R{1}(1)+rep(2)*R{2}(1),rns(2)).'; ... %"right"
    linspace(rep(1)*R{1}(1)+rep(2)*R{2}(1),rep(2)*R{2}(1),rns(1)).'; ... %"top"
    linspace(rep(2)*R{2}(1),0,rns(2)).'];                                %"left"
bnd(:,2) = [linspace(0,rep(1)*R{1}(2),rns(1)).'; ...
    linspace(rep(1)*R{1}(2),rep(1)*R{1}(2)+rep(2)*R{2}(2),rns(2)).'; ...
    linspace(rep(1)*R{1}(2)+rep(2)*R{2}(2),rep(2)*R{2}(2),rns(1)).'; ...
    linspace(rep(2)*R{2}(2),0,rns(2)).'];
bnd = bsxfun(@minus,bnd, (rep(1)*R{1}+rep(2)*R{2})/2); %Re-center


%Add wiggling to preserve positions in mesh generation (reverted later) -
%note that we only wiggle in the direction of maximum cell repetition
bndwgl = [bnd(1:rns(1),:)                      - (rep(1)>rep(2))  * bsxfun(@times,[repmat([0;1]/fac,(rns(1)-1)/2,1);0],rotR{1});
    bnd(rns(1)+1:rns(1)+rns(2),:)              - (rep(2)>=rep(1)) * bsxfun(@times,[repmat([0;1]/fac,(rns(2)-1)/2,1);0],rotR{2});
    bnd(rns(1)+rns(2)+1:2*rns(1)+rns(2),:)     + (rep(1)>rep(2))  * bsxfun(@times,[repmat([0;1]/fac,(rns(1)-1)/2,1);0],rotR{1});
    bnd(2*rns(1)+rns(2)+1:2*(rns(1)+rns(2)),:) + (rep(2)>=rep(1)) * bsxfun(@times,[repmat([0;1]/fac,(rns(2)-1)/2,1);0],rotR{2})];

bnd = unique(bnd,'rows','stable'); %Remove superfluous boundary points
bndwgl = unique(bndwgl,'rows','stable'); %Remove superfluous boundary points

%Node and connections of boundary
nodes = bndwgl;
cnct = [1:size(bndwgl,1);2:size(bndwgl,1),1].'; %Connections of unit cell boundary

%----- ADD INCLUSIONS INSIDE THE UNIT CELL -----
nodesincl = cell(rep2,1); shiftincl = ( (1:rep2) - (rep2+1)/2 );
if ~isempty(inclusion)
    for aa = 1:numel(shiftincl)
        nodesincl{aa} = bsxfun(@plus,inclusion, shiftincl(aa)*R{2});
    end

    %----- COLLISION BETWEEN BOUNDARY AND INCLUSIONS ARE DETECTED AND MERGED -----
    
    %Find all bnd nodes that are in or on the bottom (bot) or top (top) inclusion domain
    bndinincl_bot = inpolygon(bnd(:,1),bnd(:,2),nodesincl{1}(:,1),nodesincl{1}(:,2)) | ismemberf(bnd,nodesincl{1},'rows','tol',1e-10);
    bndinincl_top = inpolygon(bnd(:,1),bnd(:,2),nodesincl{end}(:,1),nodesincl{end}(:,2)) | ismemberf(bnd,nodesincl{end},'rows','tol',1e-10);
    
    %Remove vertices from the boundary and connection table that intersect with top and bottom inclusion
    kill = []; %Inserts into the nodes array to remove unwanted entries
    for bb = 1:size(bnd,1)
        if bndinincl_bot(cnct(bb,1)) || bndinincl_top(cnct(bb,1)) %Prepare to remove nodes and _lefthand_ connections that are in/on bot/top inclusion
            kill(end+1) =  cnct(bb,1);
            cnct(bb,1) = NaN;
        end
        if bndinincl_bot(cnct(bb,2)) || bndinincl_top(cnct(bb,2)) %Prepare to remove nodes and _righthand_ connections that are in/on bot/top inclusion
            kill(end+1) = cnct(bb,2);
            cnct(bb,2) = NaN;
        end
    end
    kill = unique(kill); %Remove double entries from kill
    cnct(all(isnan(cnct)')',:) = []; %Remove all connections that are completely in the forbidden region
    nodes(kill,:) = []; %Remove all boundary nodes that are in the forbidden region
    
    %Preallocate an array for converting back and forth from 'bnd' to new, modified 'nodes'
    bnd2node = [1:size(bnd,1); 1:size(bnd,1)].';  bnd2node(kill,2) = NaN;
    %Reorder connections to account for removed nodes
    if ~isempty(kill)
        for kk = fliplr(kill) %The flip is necessary to subtract in the right order
            cnct( cnct > kk ) = cnct( cnct > kk ) - 1;
            bnd2node( bnd2node(:,2) > kk ,2 ) = bnd2node( bnd2node(:,2) > kk, 2 ) - 1;
        end
    end
    %Find the nodes from top and bottom inclusion that lie inside the original boundary
    inclinbnd_bot = inpolygon(nodesincl{1}(:,1),nodesincl{1}(:,2),bnd(:,1),bnd(:,2)-1e-10);     %Small fudge factor included due
    inclinbnd_top = inpolygon(nodesincl{end}(:,1),nodesincl{end}(:,2),bnd(:,1),bnd(:,2)+1e-10); %to floating point precision
    %Remove (and sort in _global_ counterclockwise sense) the inclusion nodes that lie outside boundary
    if any(~inclinbnd_bot);  %Bottom row is sorted in local clockwise sense with pivot at -pi/2
        cntr = sum(nodesincl{1},1)/size(nodesincl{1},1); 
        nodesincl{1} = nodesincl{1}(inclinbnd_bot,:); %Remove exterior points
        phi = cart2pol(nodesincl{1}(:,1)-cntr(1), nodesincl{1}(:,2)-cntr(2));
        phi = rem(phi+pi/2,2*pi); phi(phi<0) = phi(phi<0)+2*pi; %Sort with pivot around -pi/2
        [~,Ibot]=sort(phi,'descend'); 
        nodesincl{1} = nodesincl{1}(Ibot,:);
    end
    if any(~inclinbnd_top); %Top row is sorted in local counterclockwise sense with pivot at pi/2
        cntr = sum(nodesincl{end},1)/size(nodesincl{end},1); 
        nodesincl{end} = nodesincl{end}(inclinbnd_top,:); %Remove exterior points
        phi = cart2pol(nodesincl{end}(:,1)-cntr(1), nodesincl{end}(:,2)-cntr(2));
        phi = rem(phi-pi/2,2*pi); phi(phi<0) = phi(phi<0)+2*pi; %Sort with pivot around pi/2
        [~,Itop]=sort(phi,'descend'); 
        nodesincl{end} = nodesincl{end}(Itop,:); 
    end
    
    numpincl = zeros(numel(nodesincl),1);
    for ii = 1:numel(nodesincl);  numpincl(ii) = size(nodesincl{ii},1); end %Count number of nodes in each inclusion
    
    %Find intersection points between bnd and top and bottom inclusions; these points need to be tied together in cnct
    [nanpos(:,1),nanpos(:,2)]=find(isnan(cnct)); [~,Inan]=sort(nanpos(:,1)); nanpos=nanpos(Inan,:);

    %BOTTOM INCLUSION
    if ~isempty(nanpos)
        for nn = 1:2 %Tie the first two intersections together (bottom)
            if nn == 1
                cnct(nanpos(nn,1),nanpos(nn,2)) = size(nodes,1)+1;
            elseif nn == 2
                cnct(nanpos(nn,1),nanpos(nn,2)) = size(nodes,1) + size(nodesincl{1},1);
            end
        end
    end
    %Insert the nodes and cncts that do not collide with bnd (bot)
    cnct = [cnct; size(nodes,1)+[1:numpincl(1)-1;2:numpincl(1),].'];
    if all(inclinbnd_bot); cnct = [cnct; size(nodes,1) + [numpincl(1),1]]; end %If nothing was removed we need to close the inclusion
    nodes = [nodes; nodesincl{1}];

    %MIDDLE INCLUSIONS
    for ii = 2:numel(nodesincl)-1; %These inclusions never collide, so no intersections need handling
        cnct=[cnct; size(nodes,1)+[1:numpincl(ii);2:numpincl(ii),1].'];
        nodes = [nodes; nodesincl{ii}];
    end
    
    %TOP INCLUSION
    if ~isempty(nanpos)
        for nn = 3:4 %Tie the last two intersections together (top)
            if nn == 3
                cnct(nanpos(nn,1),nanpos(nn,2)) = size(nodes,1)+1;
            elseif nn == 4
                cnct(nanpos(nn,1),nanpos(nn,2)) = size(nodes,1) + size(nodesincl{end},1);
            end
        end
    end
    %Insert the nodes and cncts that do not collide with bnd (top)
    cnct =  [cnct; size(nodes,1)+[1:numpincl(end)-1;2:numpincl(end)].'];
    if all(inclinbnd_top); cnct = [cnct; size(nodes,1) + [numpincl(end),1]]; end %If nothing was removed we need to close the inclusion
    nodes = [nodes; nodesincl{end}];
    
else  %If the inclusion array is empty, we just need to store a few arrays in their default, unmodified form
     bnd2node = [1:size(bnd,1); 1:size(bnd,1)].'; 
     bndinincl_bot = zeros(size(bnd,1),1);
     bndinincl_top = zeros(size(bnd,1),1);
end

%----- ALL COLLISIONS ASPECTS ARE NOW HANDLED; RETURN TO MESHING -----

%% SAVE BOUNDARY NODES IF REQUESTED AS OUTPUT

if nargout == 3
    boundary.nodes = nodes; 
    for bb = 1:size(bndwgl,1)         
        wglind = find(nodes(:,1) == bndwgl(bb,1) & nodes(:,2) == bndwgl(bb,2)); %Index of boundary elements in nodes
        if ~isempty(wglind);
        boundary.nodes(wglind,:) = bnd(bb,:);  %Reconstruct without wiggling and update counter
        end
    end   
    boundary.cnct = cnct; 
    boundary.Rrib = R{1}; %Ribbon-mesh is repeating along R{1}
    boundary.interior.nodes = []; boundary.interior.cnct = []; %Save all the interior nodes for seperate use
    for ii = 2:numel(nodesincl)-1; 
        boundary.interior.cnct= [boundary.interior.cnct; size(boundary.interior.nodes,1)+[1:numpincl(ii);2:numpincl(ii),1].'];
        boundary.interior.nodes=[boundary.interior.nodes; nodesincl{ii}];
        boundary.interior.cells{ii-1} = nodesincl{ii}; 
    end
    blochmesh = []; 
    origmesh = [];
    return %If 'boundary' is requested we do not perform the remaining calculation, assuming the user is only interested in the outline of the geometry
end

% %Plot the boundary connections as "movie" 
% for bb=1:size(cnct,1)
%     plot(nodes(cnct(bb,:),1),nodes(cnct(bb,:),2),'.-r'); hold on; drawnow; pause(0.01)
% end; hold off

%% ACTUAL MESHING

%Construct a default hdata
if ~exist('hdatafun','var') || isempty(hdatafun) || strcmpi(hdatafun,'const')
    hdata.fun = @(x,y) 0.825*(2*(norm(R{1},2)+norm(R{2},2)))/Nuc*ones(size(x));
else
    hdata.fun = hdatafun; clear hdatafun;
end

%Mesh the wiggled boundary with inclusions
[p,t,stats] = mesh2d(nodes,cnct,hdata,struct('plot',false));
disp(stats)

%Find boundary elements and reverse wiggling
bndind = zeros(nnz(~isnan(bnd2node(:,2))),1); bbi = 1; %Preallocation and counter
for bb = 1:size(bnd2node,1);
    if ~isnan(bnd2node(bb,2))
        bndind(bbi) = find(p(:,1) == bndwgl(bnd2node(bb,1),1) & p(:,2) == bndwgl(bnd2node(bb,1),2)); %Index of boundary elements in nodes
        p(bndind(bbi),:) = bnd(bnd2node(bb,1),:); bbi = bbi+1;  %Reconstruct without wiggling and update counter
    end
end

%Store mesh characteristics in structure
origmesh.p = p; origmesh.t = t; origmesh.area = meshArea(p,t);


%% CREATE BLOCH MESH ON A DOUGHNUT (DIRECT POINTS TO PARTNERS)

%Find boundary pairs in bndind (tie all corners to the lower left corner as well)
%connect "left" and "right" lanes
pairs_bndind = [1,rns(1)-sum(bndinincl_bot)]; %lower left and right corner
pairs_bndind = [pairs_bndind; [(rns(1)*2-sum(bndinincl_bot)-sum(bndinincl_top)+rns(2)-1):(2*(rns(1)+rns(2)-2)-sum(bndinincl_bot)-sum(bndinincl_top)) ; (rns(1)-sum(bndinincl_bot)+rns(2)-2):-1:(rns(1)-sum(bndinincl_bot)+1)].']; %left and right lanes
pairs_bndind = [pairs_bndind; [2*rns(1)-sum(bndinincl_bot)-sum(bndinincl_top)+rns(2)-2, rns(1)-sum(bndinincl_bot) + rns(2)-1] ]; %upper left and right corner

%-------------------------------------------------------------------------%
%FROM HERE ON, WE HAVE SIMPLY COPIED THE CODE FROM 'geomPeriodicInclusion'
%WITHOUT ANY CHANGE! NO ADDITIONAL ESSENTIAL CHANGES APPEAR NECESSARY.
%-------------------------------------------------------------------------%

%Find pairs in p indices
pairs = [bndind(pairs_bndind(:,1),:),bndind(pairs_bndind(:,2),:)];

%Reassign values in t-matrix
blochmesh.t = t;
blochmesh.p = p;
blochmesh.area = origmesh.area;
for pp = 1:size(pairs,1)
    swaps.ind{pp} = find(blochmesh.t==pairs(pp,2));
    swaps.orig(pp) = pp;
    blochmesh.t(blochmesh.t==pairs(pp,2)) = pairs(pp,1);  %Connect boundary-pairs in triangulation
end

%Delete points that are now redundant from bloch.p (and readjust reference points in bloch.t)
blochmesh.p(pairs(:,2),:) = []; %Delete points from bloch.p which are no longer referenced
iedge = pairs(:,1); %Indices of
deletepoints = sort(pairs(:,2),'descend');
%Adjust indices (since blochmesh.p was reordered)
for pp = 1:size(deletepoints,1)
    blochmesh.t = blochmesh.t - (blochmesh.t>deletepoints(pp));
    iedge = iedge - (iedge>deletepoints(pp));
end

%% CONSTRUCT CONVERSION TO ORIGINAL MESH FROM BLOCH MESH

%Make a new vector that stores the deleted points and their pair-values
remesh.p = [blochmesh.p;p(pairs(:,2),:)]; %All point positions
remesh.partner = iedge; %Index positions in value arrays associated with bloch.p
remesh.t = blochmesh.t;
for ss = 1:numel(swaps.ind)
    remesh.t(swaps.ind{ss}) = size(blochmesh.p,1)+swaps.orig(ss); %Readjust entries to right points
    remesh.partner(ss,2) = size(blochmesh.p,1)+swaps.orig(ss);
end
remesh.area = origmesh.area;
remesh.values = @(f) [f;f(remesh.partner(:,1))]; %Grabbing those values for some array f

%Store remesh in blochmesh so that the _real_ positions of triangle
%vertices can be accessed (necessary for evaluation of DifferentialMatrix)
blochmesh.remesh = remesh;
blochmesh.R = R;

%Store the ribbon lattice vector (there is only one!), which indicates the
%direction in which the ribbon is assumed repeating. Also store the
%associated reciprocal lattice vector Grib
if rep(1)>rep(2)                    %Ribbon-mesh is repeating along R{1}
    blochmesh.Rrib = R{2};
    blochmesh.Grib = 2*pi*R{2}; %1D: recip lat vecs are simple
elseif rep(2)>=rep(1)                 %Ribbon-mesh is repeating along R{2}
    blochmesh.Rrib = R{1};
    blochmesh.Grib = 2*pi*R{1}; %1D: recip lat vecs are simple
end

%% PLOT MESH (UNLESS REQUEST = 0)

%Default plottype
if ~exist('showmesh','var') || isempty(showmesh) || all(showmesh == 1)
    showmesh = 'remesh';
end


if any(showmesh ~= 0) %Unless requested not to plot
    cols = flatcolors(); %Colors
    
    %Figure
    set_figsize([],19,19); subaxis(1,1,1,'MT',.04,'MB',.07,'ML',.07,'MR',.04)
    switch showmesh
        case 'blochmesh'
            trimesh(blochmesh.t,blochmesh.p(:,1),blochmesh.p(:,2),'Color',cols{4});
        case 'remesh'
            trimesh(remesh.t,remesh.p(:,1),remesh.p(:,2),'Color',cols{4})
        case 'origmesh'
            trimesh(origmesh.t,origmesh.p(:,1),origmesh.p(:,2),'Color',cols{4})
        case 'pairs'
            cols = repmat(cols,1,ceil(size(pairs,1)/numel(cols))); %Make sure we have enough colors for all pairs
            for pp=1:size(pairs,1)
                plot(p(pairs(pp,:),1),p(pairs(pp,:),2),'.-','Color',cols{pp},'MarkerSize',7);hold on
            end
    end
    hold on
    for pp=1:size(pairs,1)
        plot(p(pairs(pp,:),1),p(pairs(pp,:),2),'o','color',cols{14},'Markersize',6)
    end
    %plot(bnd(:,1),bnd(:,2),'o','color',cols{14},'Markersize',6); hold off
    axis equal
    set(gca,'Fontsize',8)
    xlim([-1.075,1.075]*abs(rep(1)*R{1}(1)+rep(2)*R{2}(1))/2)
    ylim([-1.075,1.075]*abs(rep(1)*R{1}(2)+rep(2)*R{2}(2))/2)
    title(upper(showmesh),'Fontsize',9)
    drawnow;
end
