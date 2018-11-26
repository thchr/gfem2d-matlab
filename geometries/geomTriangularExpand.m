%% MESHING & MOMENTA
close all; clc; clear; cols=flatcolors;
Ns = 118; Ncirc = 70;
fprintf('\n|----- MESHING UNIT CELL -----|\n')
adivd=2.1;
theta = linspace(0,2*pi,Ncirc+1); theta = theta(1:end-1); %Inclusion boundary
hole = 1/(2*adivd)*[cos(theta);sin(theta)].';        %----||----

R = calcDirect(latticetype);
expandshrink = -.515;
inclusions=cell(1,4); 
cntrs = [-1,-1; 1,-1; 1,1; -1,1]/4;
cntrs = cntrs(:,1)*R{1} + cntrs(:,2)*R{2};
[cntrstheta,cntrsrho] = cart2pol(cntrs(:,1),cntrs(:,2));
newcntrs = bsxfun(@times, [cos(cntrstheta),sin(cntrstheta)], cntrsrho*(1+expandshrink));
for ii = 1:4
    inclusions{ii} = bsxfun(@plus, hole/2, newcntrs(ii,:));
end
mesh = geomPeriodicInclusion(latticetype,inclusions,Ns,@(x,y) .015*250/Ns*ones(size(x)),1);



%%
uc = [-R{1}-R{2};  -R{2}+R{1};  R{2}+R{1}; R{2}-R{1}]/2;

%Plot the unit cell
bzcol = cols{4}*.15+cols{8}*.85;  dotcol = 'w';
figure
patch(uc([1:end,1],1),uc([1:end,1],2),bzcol,'linewidth',.75)
hold on
for ii=1:numel(inclusions)
    patch(inclusions{ii}([1:end,1],1),inclusions{ii}([1:end,1],2),dotcol)
    plot(cntrs(:,1),cntrs(:,2),'o','color',cols{14})
    plot(newcntrs(:,1),newcntrs(:,2),'*','color',cols{21})
end; hold off
set(gca, 'box', 'off','color','w'); axis equal tight off
drawnow



