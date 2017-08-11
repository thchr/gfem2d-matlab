function vals = minmax(x)
%Compute the minimal and maximal value of an array x, possibly multivalued,
%but always real. Output is then vals = [min(x(:)),max(x(:))];

vals = [min(x(:)),max(x(:))];