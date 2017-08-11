function [blochmesh,origmesh,bnd] = geomRibbonInclusion(R,rep,inclusion,Nuc,hdatafun,showmesh)
%CALL:  [blochmesh,origmesh,bnd] =
%            geomRibbonInclusion(R,rep,inclusion,Nuc,hdatafun,showmesh)
%DESCRIPTION: Creates a folded ring-like mesh for calculations on 1D
%             Bravais lattices, as suitable for a ribbon-like restriction
%             of a 2D Bravais lattice described by lattice vectors R. The
%             functions calcRibbonCoulomb and calcLatticeDifferential
%             require as input the blochmesh structure generated herein.
%INPUT: R | Two lattice vectors as cell array of 2D row vectors [default:
%           square lattice with lattice constant 1]
%       rep | [1 x 2] vector which indicates the number of times each unit
%             cell is repeated in the ribbon-construction. Minimum value is
%             unity. rep(i) indicates repetition along R{i}. Ribbon is 
%             taken as infinite along R{i} for maximum rep(i) entry posi.
%       inclusion | [N x 2] column vector of points inside unit cell which
%                   indicate an inclusion [default: 1/4 radius disk].
%       Nuc | Total number of points along the unit cell [default: 250]
%       hdatafun | Meshing characteristics, supplied as function or string
%                  ('const' or 'antidot') [default = 'const']
%       showmesh | Set to 0 to avoid display of mesh, can otherwise be any
%                  of the strings 'blochmesh', 'remesh', 'origmesh', and
%                  'pairs'. [default = 'remesh']
%OUTPUT: blochmesh | The folded mesh supplied as a structure with fields
%                    p, t, area, remesh, and Rrib
%        origmesh | The original mesh (structure w. fields p, t, area).
%BY:                                   Thomas Christensen (22 April, 2016)

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

%Default total number of points along unit cell
if ~exist('Nuc','var') || isempty(Nuc)
    Nuc = 100;
end

%Number of points on each side from Nuc
ns = round([norm(R{1},2), norm(R{2},2)]/(2*(norm(R{1},2)+norm(R{2},2)))*Nuc);
rns = rep.*ns; %Repeated version, which accounts for the number of ribbon cells
for nn = 1:2   %Ensure that the number is odd
    if mod(rns(nn),2) == 0
        rns(nn) = rns(nn)+1;
    end
end


%Default inclusion and checking of duplicates in inclusion
if ~exist('inclusion','var') || strcmpi(inclusion,'disk') %Default inclusion
    ncirc = 35;                                          %is a 1/4 radius disk
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

%% ACTUAL MESHING

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

%Add inclusions inside the unit cell
if ~isempty(inclusion)
    for aa = ( (1:rep(1)) - (rep(1)+1)/2 )
        for bb = ( (1:rep(2)) - (rep(2)+1)/2 )
            if iscell(inclusion)                %For several disjoint inclusions
                for ii=1:numel(inclusion)
                    numpincl = size(inclusion{ii},1);
                    cnct = [cnct; +size(nodes,1)+[1:numpincl;2:numpincl,1].'];
                    nodes = [nodes; bsxfun(@plus,inclusion{ii}, aa*R{1} + bb*R{2})];
                end
            else                                %For a single inclusion
                numpincl = size(inclusion,1);
                cnct = [cnct; +size(nodes,1)+[1:numpincl;2:numpincl,1].'];
                nodes = [nodes; bsxfun(@plus,inclusion, aa*R{1} + bb*R{2})];
            end
        end
    end
end

%Construct a default hdata
if ~exist('hdatafun','var') || isempty(hdatafun) || strcmpi(hdatafun,'const')
    hdata.fun = @(x,y) 0.825*(2*(norm(R{1},2)+norm(R{2},2)))/Nuc*ones(size(x));
else
    hdata.fun = hdatafun; clear hdatafun;
end

%Mesh the wiggled boundary with inclusions
[p,t,stats] = mesh2d(nodes,cnct,hdata,struct('plot',false,'edgerefine',1));
disp(stats)

%Find boundary elements and reverse wiggling
bndind = zeros(size(bnd,1),1);
for bb = 1:size(bnd,1)
    bndind(bb) = find(p(:,1) == bndwgl(bb,1) & p(:,2) == bndwgl(bb,2)); %Position of boundary elements
    p(bndind(bb),:) = bnd(bb,:); %Reconstruct without wiggling
    if isempty(bndind(bb))
        warning('A boundary node was missed')
    end
end

%Store mesh characteristics in structure
origmesh.p = p; origmesh.t = t; origmesh.area = meshArea(p,t);


%% CREATE BLOCH MESH ON A DOUGHNUT (DIRECT POINTS TO PARTNERS)

%Find boundary pairs in bndind (tie all corners to the lower left corner as well)
if rep(1)>rep(2)                         %connect "top" and "bottom" lanes
    pairs_bndind = [1:rns(1) ; (2*rns(1)+rns(2)-2):-1:(rns(1)+rns(2)-1)].' ; %lower and upper lanes including corners
elseif rep(2)>=rep(1)                    %connect "left" and "right" lanes
    pairs_bndind = [1,rns(1)]; %lower left and right corner
    pairs_bndind = [pairs_bndind; [(rns(1)*2+rns(2)-1):(2*(rns(1)+rns(2)-2)) ; (rns(1)+rns(2)-2):-1:(rns(1)+1)].']; %left and right lanes
    pairs_bndind = [pairs_bndind; [2*rns(1)+rns(2)-2, rns(1) + rns(2)-1] ]; %upper left and right corner
end

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
