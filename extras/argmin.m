function solution=argmin(x,Constraints,Objective,options);

if isempty(x)
    x = recover(unique([depends(Constraints);depends(Objective)]));
end
if nargin < 4
    options = sdpsettings('verbose',0);
end
solvesdp(Constraints,Objective,options);
solution = double(x);
