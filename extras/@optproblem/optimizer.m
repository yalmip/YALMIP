function P = optimizer(self,x,u)
%OPTIMIZER  Container for optimization problem
%
%   OPT = OPTIMIZER(Problem,x,u) exports an object that contains
%   precompiled numerical data to be solved for varying arguments 
%   x, returning the optimal value of the expression u.
%
%   SEE OPTIMIZER for more info
%
%   Example
%
%    The following problem creates an LP with varying upper and lower
%    bounds on the decision variable.
%
%    The optimizing argument is obtained by indexing (with {}) the optimizer 
%    object with the point of interest. The argument should be a column
%    vector (if the argument has a width larger than 1, YALMIP assumes that
%    the optimal solution should be computed in several points) 
%   
%     A = randn(10,3);
%     b = rand(10,1)*19;
%     c = randn(3,1);
%
%     z = sdpvar(3,1);
%     sdpvar UB LB
%
%     Constraints = [A*z < b, LB < z < UB];
%     Objective = c'*z;
%     P = optproblem(Constraints,Objective);
%     % We want the optimal z as a function of [LB;UB]
%     optZ = optimizer(P,[LB; UB],z);
%     
%     % Compute the optimal z when LB=1, UB = 3;
%     zopt = optZ{[1; 3]}
%
%     % Compute two solutions, one for (LB,UB) [1;3] and one for (LB,UB) [2;6]
%     zopt = optZ{[[1; 3], [2;6]]}
%
%     A second output argument can be used to catch infeasibility
%     [zopt,infeasible] = optZ{[1; 3]}


P = optimizer(self.Constraints,self.Objective,self.Options,x,u);