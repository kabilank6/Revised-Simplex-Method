function [leave] = findleave(IBmatrix, as, xb, indices, n, phaseOne)

% Given entering column as and vector xb of basic variables
% findleave finds a leaving column of basis matrix Bmatrix
% It returns 0 if no column can be found (i.e. unbounded)
% leave=p indicates that the pth column of B leaves basis

% We must first calculate the most important vector in the whole world -
% this gives us the change in basic variables which are necessary to satisfy 
% Ax=b, given that we increase the nonbasic variable xs by 1. The ratio
% between xb and this vector also is important throughout - so computing it
% at the beginning is a good idea.

mostImportant = IBmatrix * as;           
ratio = xb ./ mostImportant;             

% Here we will check if we are in Phase Two and we are looking at an Artificial Variable
% as the Extended Leaving Variable Algorithm will be important for such a case. 
% We have our boolean variable to check which phase we are in and our
% artificial variables of course have indices of above n.


if((sum(((mostImportant ~= 0 )&(indices > n))) && (~phaseOne)))
    [~,leave] = max((mostImportant ~= 0)&(indices > n));
    return
end


% Otherwise we use the normal leaving variable algorithm.
% If we cannot find any values for the entering variable to be <= 0, then
% our problem must be unbounded and we must set leave to 0 accordingly.

if ~(any(mostImportant(mostImportant > 0)))     
    leave = 0;
else
    
    % We cannot choose anything with a negative ratio, or a denominator
    % that is not positive. Hence, what we will do here is set them to
    % infinity, as this guarantees that they won't end up being our
    % minimum, except in the unbounded case.
    
    ratio(ratio < 0) = Inf;       
    ratio(mostImportant <= 0) = Inf;    
    [~, leave] = min(ratio);
    
end

end