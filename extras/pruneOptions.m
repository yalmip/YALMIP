function options = pruneOptions(options)

if isequal(options.copt,options.default.copt)
    options.copt = [];
else
    [same,opsDiff] = isequalInf(options.copt,options.default.copt);
    if ~same
        options.copt = opsDiff;
    end
end
if isequal(options.xpress,options.default.xpress)
    options.xpress = [];
else
    [same,opsDiff] = isequalInf(options.xpress,options.default.xpress);
    if ~same
        options.xpress = opsDiff;
    end
end
if isequal(options.mosek,options.default.mosek)
    options.mosek = [];
else
    [same,opsDiff] = isequalInf(options.mosek,options.default.mosek);
    if ~same
        options.mosek = opsDiff;
    end
end
if isequal(options.cplex,options.default.cplex)
    options.cplex = [];
else
    [same,opsDiff] = isequalInf(options.cplex,options.default.cplex);
    if ~same
        options.cplex = opsDiff;
    end
end
if isequal(options.gurobi,options.default.gurobi)
    options.gurobi = [];
else
    [same,opsDiff] = isequalInf(options.gurobi,options.default.gurobi);
    if ~same
        options.gurobi = opsDiff;
    end
end
if isequal(options.osqp,options.default.osqp)
    options.osqp = [];
else
    [same,opsDiff] = isequalInf(options.osqp,options.default.osqp);
    if ~same
        options.osqp = opsDiff;
    end
end

function [same,d] = isequalInf(a,b)
same = 1;
if ~isempty(a)
    anames = fieldnames(a);
else
    anames = [];
end
if ~isempty(b)
    bnames = fieldnames(b);
else
    bnames = [];
end
d = [];
for i = 1:length(anames)
    e = getfield(a,anames{i});
    f = getfield(b,anames{i});
    if ~isequal(e,f)
        if isa(e,'double') && isa(f,'double') && isinf(e) && isinf(f)
            % Silly case, they're the same
        elseif isa(e,'struct') && isa(f,'struct')
            [sameHere,dnew] = isequalInf(e,f);
            if ~sameHere
                d = setfield(d,anames{i},dnew);
                same = same & sameHere;
            end
        else
            same = 0;
            d = setfield(d,anames{i},e);
        end
    end
end


