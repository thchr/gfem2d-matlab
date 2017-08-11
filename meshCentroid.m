function centroid = meshCentroid(p,t)
%CALL:        centroid = meshCentroid(p,t)
%DESCRIPTION: Calculation the centroids of triangles in a mesh (p,t)

centroid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:) ) / 3;


