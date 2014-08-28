function [feasible,feaslistLMI] = isfeasible(F,tol)

if nargin == 1
    tol = 0;
end
feaslistLMI = checkset(F);
feasible = all(feaslistLMI >= -tol);
