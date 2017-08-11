clear all; close all; clc;

[cols,nams] = flatcolors;
R{1} = [1,0]; R{2} = [0,1];
nR{1} = R{1}(:)/norm(R{1},2); nR{2} = R{2}(:)/norm(R{2},2);

nside = 51;
r = 1;
corners = [R{1} + R{2}; R{1} - R{2}; -R{1} - R{2}; -R{1} + R{2}]/2;
fac=100;
%Outer boundary (square)
node(:,1) = [linspace(0,R{1}(1),nside).'; ...
             linspace(R{1}(1),R{1}(1)+R{2}(1),nside).' + [repmat([0;1]/fac,(nside-1)/2,1);0]; ...
             linspace(R{1}(1)+R{2}(1),R{2}(1),nside).'; ...
             linspace(R{2}(1),0,nside).' - [repmat([0;1]/fac,(nside-1)/2,1);0]];
node(:,2) = [linspace(0,R{1}(2),nside).' - [repmat([0;1]/fac,(nside-1)/2,1);0]; ...
             linspace(R{1}(2),R{1}(2)+R{2}(2),nside).'; ...
             linspace(R{1}(2)+R{2}(2),R{2}(2),nside).' + [repmat([0;1]/fac,(nside-1)/2,1);0]; ...
             linspace(R{2}(2),0,nside).'];
node = unique(node,'rows','stable'); %Remove superfluous boundary points
node = node/2;

cnct = [1:size(node,1);2:size(node,1),1].';
% Make mesh
[ptemp,t] = mesh2d(node,cnct); 
%break
%Extrapolate mesh to the remaining regions
M = [nR{1}(1) , nR{2}(1); 
    nR{1}(2) , nR{2}(2)]; invM = inv(M);
pnat = (M\ptemp.').';

pNE = ptemp; tNE = t; 
pSE = [pnat(:,1)*nR{1}(1)-pnat(:,2)*nR{2}(1), pnat(:,1)*nR{1}(2)-pnat(:,2)*nR{2}(2)];
pNW = [-pnat(:,1)*nR{1}(1)+pnat(:,2)*nR{2}(1), -pnat(:,1)*nR{1}(2)+pnat(:,2)*nR{2}(2)];
pSW = [-pnat(:,1)*nR{1}(1)-pnat(:,2)*nR{2}(1), -pnat(:,1)*nR{1}(2)-pnat(:,2)*nR{2}(2)];

%Create a new mesh without duplicate nodes
p = [pNE; pSE; pNW; pSW]; 
t = [t;t+size(ptemp,1);t+2*size(ptemp,1);t+3*size(ptemp,1)];
[p,t] = fixmesh(p,t); %Remove duplicates
[~,~,~,bnd] = connectivity(p,t);
figure
trimesh(t,p(:,1),p(:,2),'color',[0,0,0]); axis equal; hold on;
%plot(p(bnd,1),p(bnd,2),'s'); hold off;
ie = find(bnd); pe = p(ie,:);
plot(pe(:,1),pe(:,2),'s'); 
hold off

%%
tol = 1e-8;
pmR1 = bsxfun(@minus,pe, R{1});
X=abs(bsxfun(@minus,pmR1(:,1),pe(:,1).'))< tol;
Y=abs(bsxfun(@minus,pmR1(:,2),pe(:,2).'))< tol;

[top1,bot1] = find(X&Y);

pmR2 = bsxfun(@minus,pe, R{2});
X=abs(bsxfun(@minus,pmR2(:,1),pe(:,1).'))< tol;
Y=abs(bsxfun(@minus,pmR2(:,2),pe(:,2).'))< tol;
[top2,bot2] = find(X&Y);

figure
plot(pe(top1,1),pe(top1,2),'.r'); hold on
plot(pe(bot1,1),pe(bot1,2),'.m');
plot(pe(top2,1),pe(top2,2),'sg');
plot(pe(bot2,1),pe(bot2,2),'s');
%% FIGURE
set_figsize([],14,14)
hold on;
%Unit cell
plot(corners(:,1),corners(:,2),'o','MarkerSize',6.5,'MarkerFaceColor',[.2,.2,.2],'MarkerEdgeColor',[.7,.7,.7]);
plot(corners([1:end,1],1),corners([1:end,1],2),':','Color',[.7,.7,.7])
%Boundary nodes of original square
plot(node(:,1),node(:,2),'.','color',cols{8},'MarkerSize',15)
%Lattice vectors
plot([0,R{1}(1)],[0,R{1}(2)],'Color',cols{14});
plot([0,R{2}(1)],[0,R{2}(2)],'Color',cols{14})
%Mesh vertex points
plot(pNE(:,1),pNE(:,2),'.','color',cols{1},'MarkerSize',5)
plot(pSE(:,1),pSE(:,2),'s','color',cols{2},'MarkerSize',5)
plot(pSW(:,1),pSW(:,2),'o','color',cols{3},'MarkerSize',5)
plot(pNW(:,1),pNW(:,2),'d','color',cols{4},'MarkerSize',5)

axis equal off