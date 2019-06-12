function newConstraint =  sampledScenarioChanceFilter(b,c,distribution,One_minus_confidencelevel,w,options);

gamma = 1-One_minus_confidencelevel;
if isempty(options.chance.scenario.delta)
    N = options.chance.N;    
else
    delta = options.chance.scenario.delta;
    n = max(length(depends(b)),length(depends(b))); % FIX
    N = ceil((2/gamma)*(2*log(2/gamma) + log(1/delta))+2*n);
end
W = [];for i = 1:N;W = [W dataSampler(distribution,size(w))];end
newConstraint = [-b-c'*W <= 0];