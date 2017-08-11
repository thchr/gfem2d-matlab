function phi = calcDipoleSource(r,r0,d)
%CALL:      phi = calcDipoleSource(r,r0,d)
%DESCRIPTION: Calculates the potential, phi, due to a dipole source at
%position r0 with dipole moment d (the latter two being 3-vectors in xyz):
%       phi = [ (r-r0).d / |r-r0|^3 ]
%The dipole moment d indicates the direction and magnitude of the dipole.
%The expression can be derived by considering two charges of opposite sign,
%at positions r0 + delta and r0 - delta, respectively, with delta propto d.
%INPUT:     r  | Positions to evaluate the potential [N x 3 array]
%           r0 | Position of dipole [1 x 3 array]
%           d  | Dipole moment as vector [1 x 3 array]
%OUTPUT:   phi | The electric potential associated with a dipole d at r0,
%                evaluated at positions r. Note that units are not
%                included, and, as a consequence, the magnitude of d is
%                meaningless apart from relative magnitude.
%Thomas Christensen                                         02 May, 2016

if size(r,2) == 2; %If r is a 2D array, we assume the z-coordinate is zero
    r(:,3) = zeros(size(r,1),1);
end
R = bsxfun(@minus,r,r0); %r-r0 as [N x 3 array]
phi = R*d(:)./( R(:,1).^2 + R(:,2).^2 + R(:,3).^2 ).^(3/2);
