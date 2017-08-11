function [gradf,rc] = calcGradient(p,t,f)
%Calculates the gradient of a vertex (p) specified function (f) on a mesh
%(t), as applicable to the faces of each mesh-element.

gradf = zeros(size(t,1),2);
area = meshArea(p,t);

%ptx = p(t,1); pty = p(t,2); 
%rjk(:,1:2) = p(t(:,1),:) - p(t(:,2),:); %rab
%rjk(:,3:4) = p(t(:,2),:) - p(t(:,3),:); %rbc
%rjk(:,5:6) = p(t(:,3),:) - p(t(:,1),:); %rca
for jj = 1:size(t,1)
    gradf(jj,:) = ( f(t(jj,1))*( p(t(jj,2),:) - p(t(jj,3),:) ) + ...
                    f(t(jj,2))*( p(t(jj,3),:) - p(t(jj,1),:) ) + ...
                    f(t(jj,3))*( p(t(jj,1),:) - p(t(jj,2),:) ) )/(2*area(jj));   
end
gradf = ( [0,1;-1,0]*gradf.' ).'; %Flip items

if nargout == 2; %Calculate centroid positions to associate with gradf
    rc = meshCentroid(p,t);
end