function cpnt = calcQuadraturePoints(rl,rm,rn)
%Computes the centerpoints of 16 subtriangle of a single triangle as 
%specified by vertices rl, rm, rn. The is needed for the 16-point
%quadrature scheme employed in calcCoulombNumQuad and calcLatticeCoulomb.

%Subtriangles corners/vertex points
pnt(1, :)  = (3 * rl      + rm) / 4;
pnt(2, :)  = (3 * rl      + rn) / 4;
pnt(3, :)  = (rl          + rm) / 2;
pnt(4, :)  = (2 * rl + rm + rn) / 4;
pnt(5, :)  = (rl          + rn) / 2;
pnt(6, :)  = (rl      + 3 * rm) / 4;
pnt(7, :)  = (rl + 2 * rm + rn) / 4;
pnt(8, :)  = (rl + rm + 2 * rn) / 4;
pnt(9, :)  = (rl      + 3 * rn) / 4;
pnt(10, :) = (3 * rm      + rn) / 4;
pnt(11, :) = (rm          + rn) / 2;
pnt(12, :) = (rm      + 3 * rn) / 4;

%Center points in subtriangles
cpnt(1,:) = (rl        + pnt(1, :)  + pnt(2, :))  / 3;
cpnt(2,:) = (pnt(1, :) + pnt(3, :)  + pnt(4, :))  / 3;
cpnt(3,:) = (pnt(1, :) + pnt(4, :)  + pnt(2, :))  / 3;
cpnt(4,:) = (pnt(2, :) + pnt(4, :)  + pnt(5, :))  / 3;
cpnt(5,:) = (pnt(3, :) + pnt(6, :)  + pnt(7, :))  / 3;
cpnt(6,:) = (pnt(3, :) + pnt(7, :)  + pnt(4, :))  / 3;
cpnt(7,:) = (pnt(4, :) + pnt(7, :)  + pnt(8, :))  / 3;
cpnt(8,:) = (pnt(4, :) + pnt(8, :)  + pnt(5, :))  / 3;
cpnt(9,:) = (pnt(5, :) + pnt(8, :)  + pnt(9, :))  / 3;
cpnt(10,:) = (pnt(6, :) + rm         + pnt(10, :)) / 3;
cpnt(11,:) = (pnt(6, :) + pnt(10, :) + pnt(7, :))  / 3;
cpnt(12,:) = (pnt(7, :) + pnt(10, :) + pnt(11, :)) / 3;
cpnt(13,:) = (pnt(7, :) + pnt(11, :) + pnt(8, :))  / 3;
cpnt(14,:) = (pnt(8, :) + pnt(11, :) + pnt(12, :)) / 3;
cpnt(15,:) = (pnt(8, :) + pnt(12, :) + pnt(9, :))  / 3;
cpnt(16,:) = (pnt(9, :) + pnt(12, :) + rn)         / 3;