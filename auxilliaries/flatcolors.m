function [col,colnam] = flatcolors(sortchoice)
%Gives a cell array of pretty flat colors for use in plots

col{1}  = [ 217, 30, 24 ];    colnam{1}  = 'Thunderbird, red';
col{2}  = [ 75, 119, 190 ];   colnam{2}  = 'Steel blue';
col{3}  = [ 30, 130, 76 ];    colnam{3}  = 'Salem, green';
col{4}  = [ 20, 20, 22 ];     colnam{4}  = 'Nearly black';
col{5} = [ 242, 120, 75 ];    colnam{5}  = 'Crusta, off-orange';
col{6}  = [ 154, 18, 179 ];   colnam{6}  = 'Seance, purple';
col{7}  = [ 247, 202, 24 ];   colnam{7}  = 'Ripe lemon, yellow';
col{8}  = [ 171, 183, 183 ];  colnam{8}  = 'Edward, light grey';
col{9} = [ 78, 205, 196 ];    colnam{9}  = 'Medium turquoise';
col{10} = [ 103, 65, 114 ];   colnam{10} = 'Light purple';
col{11} = [ 1, 152, 117 ];    colnam{11} = 'Green haze';
col{12} = [ 200, 247, 197 ];  colnam{12} = 'Madang';
col{13} = [ 189, 195, 199 ];  colnam{13} = 'Silver';
col{14} = [ 31, 58, 147 ];    colnam{14} = 'Jackson purple';
col{15} = [ 144, 198, 149 ];  colnam{15} = 'Dark green sea';
col{16} = [ 248, 148, 6 ];    colnam{16} = 'California';
col{17} = [ 101, 198, 187 ];  colnam{17} = 'Downy';
col{18} = [ 38, 166, 91 ];    colnam{18} = 'Eucalyptus';
col{19}  = [ 34, 49, 63 ];    colnam{19}  = 'Ebony clay';
col{20}  = [ 211, 84, 0 ];    colnam{20}  = 'Burnt orange';
col{21}  = [ 150, 40, 27 ];   colnam{21}  = 'Old brick';
col{22}  = [ 207, 0, 15 ];    colnam{22}  = 'Bright red';
col{23}  = [ 192, 57, 43 ];   colnam{23}  = 'Dull red';
col{24}  = [ 68,108,179 ];    colnam{24}  = 'San Marino';


col = cellfun(@(x) x/255, col,'UniformOutput',0); %Normalize to unity elements

%Sort colors according to requested scheme (if requested)
if exist('sortchoice','var') && ~isempty(sortchoice)
    if strcmpi(sortchoice,'hue')
        sortval = cellfun(@(col) atan2( sqrt(3) * (col(2)-col(3)) , 2 * col(1) - col(2) - col(3) ), col);
    elseif strcmpi(sortchoice,'brightness')
        sortval = cellfun(@(col) sqrt( sum([0.299,.587,.114].*col.^2) ), col);
    end
    [~,sortorder] = sort(sortval,'ascend');
    col = col(sortorder); colnam = colnam(sortorder);
end
