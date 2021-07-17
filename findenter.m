function [as, cs, s] = findenter(Amatrix, pi, c, indices, n)

% Given the complete m by n matrix Amatrix,  
%       the complete cost vector c with n components
%       the vector pi with m components
% findenter finds the index of the entering variable and its column 
% It returns the column as, its cost coefficient cs, and its column index s
% Returns s=0 if no entering variable can be found (i.e. optimal)
% This will happen when minimum reduced cost > tolerance
% where tolerance = -1.0E-6 

tolerance = -1.0e-6;

for i = 1 : n
    
    % We only want to perform pricing on columns that are non-basic. To do
    % this, we must only consider the columns where there are no non-zero
    % values.
    
    if ~(any((indices == i) ~= 0))
        rs = c(i, 1) - pi' * Amatrix(:, i);
        if (rs < tolerance)
            s = i;
            as = Amatrix(:, i);
            cs = c(i, 1);
            return;
        end
    end
    
    as = 0;
    cs = 0;
    s = 0;
end

end