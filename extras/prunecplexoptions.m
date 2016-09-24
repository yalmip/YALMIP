function options = prunecplexoptions(options)

try
    if strfind(version,'2016')
        o1 = cplexoptimset;
    else
        o1 = cplexoptimset('cplex');
    end
    o2 = options.cplex;
    n = fieldnames(o1);
    for i = 1:length(n)
        if isequal(o1.(n{i}),o2.(n{i}))
            options.cplex = rmfield(options.cplex,n{i});
        end
    end
catch
end