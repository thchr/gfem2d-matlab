function [area,areavec] = meshArea(p,t)
%CALL:      area = meshArea(p,t)
%Calculates the area of triangles specified by p and t

area = triarea(p,t)/2; %NOTE THE FACTOR OF 2!

if nargout == 2 %This is a (row) vector which, when dotted onto a vertex-specified function, integrates that function over (p,t)
    areavec = zeros(1,size(p,1));
    for jj = 1:size(t,1)
        areavec(t(jj,:)) = areavec(t(jj,:)) + area(jj)/3;
    end
end