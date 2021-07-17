function [z, x, pi, indices, exitflag] = rsm(A, b, c, m, n)

% This function solves a linear program in standard computational form
% (that is, in the form of min cx s.t. Ax=b, x>=0) using the revised
% simplex method. The exitflag will be 0 if we successfully solved the
% problem, 1 if infeasible, and -1 if unbounded.

% We need to start off by initialising our variables. We start with the
% phase one simplex method. We can use the same function for phase one and
% phase two with little differences (except for how we find the leaving
% variable). Ultimately, the way for us to check whether we are in Phase
% One or Two will be to check if any of our 'real' or non-artificial
% variables have a non-zero value within the simplex method. If this is the
% case, we are in phase two, as all real variables must be zero in phase
% one.

oneArray = ones(m, 1);
nZeroArray = zeros(n, 1);
costOne = [nZeroArray; oneArray];
IBmatrix = eye(m);
indices = [n + 1: 1: n + m]';



[z, x,pi, indices, ~] = simplex(A, b, costOne, m, n, IBmatrix, indices);

% Our cost is only feasible if z == 0 in phase One as all 'real' variables
% must have a cost of 0 and all artificial variables should have a value of
% 0 (that is, at optimality z* == 0). If this is not the case, we know that our 
% problem is infeasible (exitflag = 1) and we
% can terminate the algorithm as Phase Two is not necessary.

if (z ~= 0)
    exitflag = 1;
    return;
end

% We have not left the algorithm at this stage. This means that our problem
% is at optimality. As such, we can proceed with the Phase Two Simplex
% method. 

mZeroArray = zeros(m, 1);
A2 = [A, IBmatrix];
phaseTwoInverseBasis = inv(A2(:, indices));
costTwo = [c; mZeroArray];

[z, x, pi, indices, exitflag] = simplex(A, b, costTwo, m, n, phaseTwoInverseBasis, indices);

end