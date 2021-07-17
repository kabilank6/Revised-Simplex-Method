function [z, x, pi, indices, exitflag] = simplex(A, b, c, m, n, IBmatrix, indices)

% Solves min cx s.t. Ax=b, x>=0
% starting with basic variables listed in vector indices
% and basis matrix Bmatrix
% exitflag is 1 if solved successfully and -1 if unbounded
% returns optimal vector x and its value z, along with pi, and indices of basic variables


xOne = zeros(m + n, 1);
cb = c(indices, 1);
xb = IBmatrix * b;
basicCosts = c(1 : n, 1);


while true
    
    % Our pi vector can help with maintaining excellent computational
    % efficiency upon finding the entering variable. As such, we will
    % compute it here.
    
    pi = (cb' * IBmatrix)';
    
    % We will find our entering variable using findenter.
    
    [as,cs,s] = findenter(A, pi, basicCosts, indices, n);
    
    % We have a phaseOne variable to track whether we are in phase one. If
    % this is false we are in phase two. We initialise it to true as we
    % perform a check later to see whether it is false. 
    
    phaseOne = true;
    if ~ s
        xOne(indices) = xb;
        x = xOne(1 : n, 1);
        z = c' * xOne;
        exitflag = 0;       
        break
    end
    
    
    % If we are in Phase One, then all of our costs for real
    % (non-artificial) variables will be zero. If any of these costs are not
    % equal to zero then we are in Phase Two. As such, we need to update
    % our boolean variable accordingly.
    
    if any(c(1: n, 1)) 
        phaseOne = false;
    end
    
    % Next step is finding the leaving variable. It does matter whether we
    % are in Phase One or Phase Two here, so our boolean comes in handy. 
    
    [leave] = findleave(IBmatrix, as, xb, indices, n, phaseOne);
    
    % If leave == 0, we don't have a leaving variable. Subsequently, we can
    % conclude that our program is unbounded. As such, we must set exitflag
    % to -1.
    
    if ~ leave
        xOne(indices) = xb;
        x = xOne(1:n,1);
        z = c' * xOne;
        exitflag = -1;      
        break
    end
    
    % Otherwise, we can perform the Gauss-Jordan pivot.
    
    [IBmatrix, indices, cb,xb] = updateGJ(IBmatrix, indices, cb, cs, as, s, leave, xb);
end

end