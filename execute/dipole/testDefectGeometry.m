clear all; clc; close all;
addpath(genpath('../../'));

ConstantsUnits0; [cols,nams]=flatcolors();

cut = .4;
repincl =  [11,7]; %[13,9];% ;
Ncirc = 29;  %28
adivd = 2;
bndshape = 'rectangle';
auxmesh = 'zoom';
usedefect = 'yes';
collisiontest = 0;
R{1} = [1,0]; R{2} = [cosd(60),sind(60)];

rep = repincl-cut*2;
repincl = bsxfun(@minus,repincl*2,mod(repincl,2));

%% ----- SETUP OF BOUNDARY AND MESHING -----
%boundary of the ribbon unit cell, which consist of several lattice unit cells (first x-coords, then y-coords)
switch bndshape
    case 'trapezoid'
        bnd = [zeros(1,2); rep(2)*R{2}; rep(1)*R{1}+rep(2)*R{2}; rep(1)*R{1}];
    case 'rectangle'
        bnd = [zeros(1,2); [0,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)]; rep(1)*R{1} + rep(2)*R{2}; [rep(1)*R{1}(:,1) + rep(2)*R{2}(:,1),0]];
    case 'defectrect'
        w=1-cut/2; h = 2*R{2}(2); repcut = ceil(rep(1)*.75)+(1-cut)/adivd/2 + 0;
        bnd = [zeros(1,2);
            [0,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)];
            [repcut*R{1}(:,1) - w/2,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)]; [repcut*R{1}(:,1) - w/2,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2) - h/50];
            [repcut*R{1}(:,1) - w/2,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2) - h];
            [repcut*R{1}(:,1) + w/2,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2) - h]; [repcut*R{1}(:,1) + w/2,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2) - h/50]
            [repcut*R{1}(:,1) + w/2,rep(1)*R{1}(:,2) + rep(2)*R{2}(:,2)];
            rep(1)*R{1} + rep(2)*R{2};
            [rep(1)*R{1}(:,1) + rep(2)*R{2}(:,1),0]];
end
bnd = bsxfun(@minus,bnd, (rep(1)*R{1}+rep(2)*R{2})/2); %Re-center
bnd = unique(bnd,'rows','stable'); %Remove superfluous boundary points

%inclusion basic shape
theta = linspace(0,2*pi,Ncirc); theta=theta(1:end-1).'+theta(2)/2;
circ = [cos(theta),sin(theta)]/(2*adivd);

%add inclusions inside the unit cell
ii = 1;
for aa = ( (1:repincl(1)) - (repincl(1)+1)/2 )
    for bb = ( (1:repincl(2)) - (repincl(2)+1)/2 )
        incl{ii} = bsxfun(@plus,circ, aa*R{1} + bb*R{2});
        ii = ii + 1;
    end
end

%make a defect (this entire thing is prone to be slightly buggy)
if strcmpi(usedefect,'yes') || all(usedefect == 1)
    w = 1; h = 2; Nd = 3;
    defect(:,1) = [linspace(w/2,-w/2,2).'; linspace(-w/2,-w/2,Nd).';  linspace(-w/2,w/2,2).';   linspace(w/2,w/2,Nd).' ];
    defect(:,2) = [linspace(h,h,2).'; linspace(h,-h,Nd).'; linspace(-h,-h,2).'; linspace(-h,h,Nd).'];
    defect = unique(defect,'rows','stable');
    defect = bsxfun(@plus,0.5*R{1} + [0,rep(2)*R{2}(:,2)/2],defect); %Move to boundary
    
    origdefect = defect;
    
    %test if any inclusions overlap with defect; if so, we have to handle that
    count = 1;
    for ii = 1:numel(incl)
        indefect = inpolygon(incl{ii}(:,1),incl{ii}(:,2),defect(:,1),defect(:,2));
        if any(indefect);
            rmvincl(count).cellind = ii;
            rmvincl(count).arrayind = find(~indefect);
            count = count + 1;
        end
    end
    for rr = 1:numel(rmvincl); %Remove points from defect list which lie within an inclusion
        definincl = inpolygon(defect(:,1),defect(:,2),incl{rmvincl(rr).cellind}(:,1),incl{rmvincl(rr).cellind}(:,2));
        defect = defect(~definincl,:);
    end
    if collisiontest == 1 %This entire thing is buggy under many circumstances, especially when it intersects in "non-convex hull"-like ways
        for rr = 1:numel(rmvincl) %Add the inclusions to the defect list
            defect = [defect; incl{rmvincl(rr).cellind}(rmvincl(rr).arrayind,:)];
        end
        
        unsort = defect(2:end,:); defect = defect(1,:);
        for dd = 1:size(unsort,1) %Sort the new defect list
            temp = unsort; conflict = 1;
            while conflict == 1 && ~all(isnan(temp(:))) %We don't want lines that go across the original defect
                [~,I] = min(pdist2(defect(dd,:),temp));
                midpoint = (defect(dd,:)+temp(I,:))/2;
                [IN,ON] = inpolygon(midpoint(:,1),midpoint(:,2),origdefect(:,1),origdefect(:,2));
                conflict = IN*~ON;
                if conflict == 1
                    for rr = 1:numel(rmvincl)
                        midpoint = (defect(dd,:)+temp(I,:))/2;
                        [IN,ON] = inpolygon(midpoint(:,1),midpoint(:,2),incl{rmvincl(rr).cellind}(:,1),incl{rmvincl(rr).cellind}(:,2));
                        conflict = IN*~ON;
                        if conflict == 0; break; end
                    end
                end
                hold off; axis equal; drawnow; pause
                temp(I,:) = NaN;
            end
            defect(dd+1,:) = unsort(I,:);
            unsort(I,:) = [];
        end
    end
    %Remove the inclusions that overlapped with the original defect
    for rr = numel(rmvincl):-1:1;  incl(rmvincl(rr).cellind) = [];  end
    
    %Finally, add the new composite defect to the inclusion list
    incl{end+1} = defect;
    
    %Text for saving
    defecttxt = 'WithDefect';
else
    defecttxt = 'NoDefect';
end




%mesh the structure
[p,t,nodes,cnct] = geomFancyInclusion(bnd,incl,@(x,y) 150*ones(size(x)),'mesh');