function options = pruneOptions(options)

% Mosek is easy to check since there are no inf parameters...
if isequal(options.mosek,options.default.mosek)
    options.mosek = [];
end
if isequal(options.cplex,options.default.cplex)
    options.cplex = [];
end
if isequal(options.gurobi,options.default.gurobi)
    options.gurobi = [];
end


function [same,d] = isequalInf(a,b)
same = 1;
anames = fieldnames(a);
bnames = fieldnames(b);
for i = 1:length(anames)
    e = getfield(a,anames{i});
    f = getfield(b,anames{i});
    if ~isequal(e,f)
        if isinf(e) && isinf(f)
            
        elseif isa(a,'struct') && isa(b,'struct')
            sameHere = isequalInf(e,f);
            same = same & sameHere;
        else
            same = 0;
            setfield(d,anames{i},b);
        end
    end
end


