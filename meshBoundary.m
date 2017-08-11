function [Ip,pedge,It,tedge,tedge2] = findBoundary(p,t)
%CALL:          [I,pedge,It,tedge,tedge2] = findBoundary(p,t)
%DESCRIPTION: Calculates the indices of the edge vertices for a triangular
%mesh (p,t) using an adjancency matrix approach. Outputs the corresponding
%indices Ip (intended to be indexed into p) and the boundary points pedge.
%In addition, the face-indices It (to be indexed into t) that correspond to
%triangles on the boundary, as well as an associated vector tedge = t(It,:)
%METHOD: Copied from the blog https://mathproblems123.wordpress.com/ ...
%2015/04/21/%identifying-edges-and-boundary-points-2d-mesh-matlab/ - the
%algorithm works by identifying edges which are sides to only one triangle
%                                                      TC | 20 April, 2016

% number of vertices in mesh
np = size(p,1);
% adjacency matrix
A = min(sparse(t(:,1),t(:,2),1,np,np)+sparse(t(:,2),t(:,3),1,np,np)+sparse(t(:,3),t(:,1),1,np,np),1);
A = min(A+A',1);
% this finds the boundary points, whose indices are stored in I
B = A^2.*A==1;
Ip = find(sum(B,2)>0);

if nargout >= 2 % calculate the boundary points
    pedge = p(Ip,:);
end

if nargout >= 3 % my own (ugly) addition to the above; finds the face-indices into t
    It = zeros(numel(Ip)-1,1); count = 1;
    for tt=1:size(t,1);
        ch = any(t(tt,1) == Ip) + any(t(tt,2) == Ip) + any(t(tt,3) == Ip);
        if ch == 2 %If there are two nodes on the edge, it must be a proper edge-face
            It(count) = tt;
            count = count+1;
        end
    end
end

if nargout >= 4 % find the face-indices into p (i.e. a new "edge-only" t)
    tedge = t(It,:);
end

if nargout >= 5 % similar to tedge, but with edges vertices only (Nx2)
    tedge2 = zeros(size(tedge,1),2); 
    for ee=1:size(tedge,1)
        count2 = 1;
        for vv = 1:3
            if any(tedge(ee,vv) == Ip)
                tedge2(ee,count2) = tedge(ee,vv); count2 = count2+1;
            end
        end
    end
end
