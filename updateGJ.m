function [IBmatrix, indices, cb,xb] = updateGJ(IBmatrix, indices, cb, cs, as, s, leave,xb)


mostImportant = IBmatrix * as;
GJ = [xb, IBmatrix, mostImportant];     
ratioNumerator = GJ(:, end);
ratioDenominator = GJ(leave, end);
ratio = ratioNumerator ./ ratioDenominator;   
totalLength = length(ratio);
row = GJ(leave,:);                  

% We now have the elements to conduct a Gauss-Jordan elimination.

for i = 1:totalLength
    if i ~= leave
        GJ(i, :) = GJ(i, :)- ratio(i, 1) .* row;     
    else
        GJ(i, :) = GJ(i, :) * (1 / GJ(i, end));        
    end
end

xb = GJ(:, 1);
IBmatrix = GJ(:, 2: end - 1);
indices(leave, 1) = s;
cb(leave, 1) = cs;

end