function plotVertexData(p,t,v,morecommand)
%CALL:              plotVertexData(p,t,v,morecommand)
%Plots the values v defined on the vertex points p of the triangle t by
%means of the trisurf command. The 'morecommand' value can be specified as
%a string (or a cell of strings), which will then evaluate.


set_figsize([],20,20);
trisurf(t,p(:,1),p(:,2),v,'EdgeColor','none'); 
shading interp; 
colormap('bluewhitered_mod')
view(2); 
axis equal off; 
drawnow; 

%Evaluate additional commands provided in 'morecommand'
if nargin == 4
    if ischar(morecommand)
        eval(morecommand)
    elseif iscell(morecommand)
        for cc = 1:numel(morecommand)
            eval(morecommand{cc})
        end
    end
end
