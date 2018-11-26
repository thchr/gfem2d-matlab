function R = calcDirect(lattice)
%CALL:      R = calcDirect(lattice)
%DESCRIPTION: Calculates the direct lattice vectors of a state described by
% 'lattice': returns answer as a two-element cell array 'R'

%Default lattice vectors (if R is empty): a square of side 1
if ~exist('lattice','var') || isempty(lattice) || all(strcmpi(lattice,'square'))
    R{1} = [1,0]; R{2} = [0,1];
elseif strcmpi(lattice,'triangular')
    R{1} = [1,0]; R{2} = [cosd(60),sind(60)];
elseif all(size(lattice) == [1,1]) %If lattice is a scalar, we interpret it as the
    Rl = lattice;                  %width/height of a square unit cell
    R{1} = Rl*[1,0]; R{2} = Rl*[0,1];
end