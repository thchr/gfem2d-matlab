function intf = integrateMeshFunction(p,t,f,areavec)
%CALL:          fint = integrateMeshFunction(p,t,f,areavec);
%Analytically integrate a vertex-specified function 'f' (expected to be in
%column form [N x 1]; if supplied as matrix [N x M], the M columns will be
%considered as distinct functions to be integrated) triangular mesh, with
%triangulation matrix 't' and vertex-coordinates 'p'. A fourth input
%'areavec' can be supplied to avoid calculating the mesh elements' areas
%multiple times

sizef = size(f); sizet = size(t);
if sizef(1) ~= sizet(1) && sizef(2) == sizet(1) %Make sure the integrand is in the right form
    f = f.';
end

%Actual analytical summation/calculation based on areavec-technique
if nargin ~= 4
    [~,areavec] = meshArea(p,t);
end
%All dot products between columns of f and areavec as vectorized operation
intf = sum(bsxfun(@times,f,areavec.'),1).';


%% OLD IMPLEMENTATION
% function intf = integrateMeshFunction(f,t,areas)
% %CALL:          fint = integrateMeshFunction(f,t,areas);
% %Analytically integrate a vertex-specified function 'f' (expected to be in
% %column form [N x 1]; if supplied as matrix [N x M], the M columns will be
% %considered as distinct functions to be integrated) triangular mesh, with
% %triangulation matrix 't' and mesh-areas 'areas'.
% 
% sizef = size(f); sizet = size(t);
% if sizef(1) ~= sizet(1) && sizef(2) == sizet(1) %Make sure the integrand is in the right form
%     f = f.';
% end
% 
% vertInt = [2,1,1;1,2,1;1,1,2]/12; %Integration matrix for vertex-specified functions
% 
% intf = zeros(size(f,2),1);
% for cc = 1:size(f,2)
%     integrand_vertices = f(:,cc);
%     integrand_trimesh = bsxfun(@times, areas,integrand_vertices(t)); 
%     intf(cc) = sum(sum(vertInt*integrand_trimesh.',1),2);
% end
