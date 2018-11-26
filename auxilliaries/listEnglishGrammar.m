function str = listEnglishGrammar(items)
%Just a simple helper function to turn e.g. [1,4,7,8] into '1, 4, 7, & 8'
%or [1,7] into '1 & 7'.

if numel(items) == 1
    str = num2str(items);
elseif numel(items) == 2
    str = [num2str(items(1)) ' & ' num2str(items(2))];
else
    str = '';
    for nn = 1:numel(items)-1
        str = [str num2str(items(nn)) ', '];
    end
    str = [str '& ' num2str(items(end))];
end
