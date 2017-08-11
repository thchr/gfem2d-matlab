function [G,A_uc] = calcReciprocal(R)
%CALL:          [G,A_uc] = calcReciprocal(R)
%Calculates the reciprocal lattice vectors G and the area of the unit cell
%spanned by two 2D lattice vectors R{1} and R{2}.
%The following variable form is expected:
%       R  : 2-element cell array with entries R{i} : [2x1] column vector
%Output:
%       G  : 2-element cell array with entries G{2} : [2x1] column vector

for ii = 1:2; %If R is supplied in row-form, we change it to column form.
    if all(size(R{ii}) == [1,2])
        R{ii} = R{ii}.';
    end
end

A_uc = norm(cross([R{1};0],[R{2};0]),2);
G{1} = (2*pi/A_uc^2) * ( cross([R{2};0],cross([R{1};0].',[R{2};0])) );
G{2} = (2*pi/A_uc^2) * ( cross([R{1};0],cross([R{2};0].',[R{1};0])) );

%Remove superfluous third element
G{1} = G{1}(1:2); G{2} = G{2}(1:2); 