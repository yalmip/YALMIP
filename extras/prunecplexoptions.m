function options = prunecplexoptions(options)

try
    o1 = cplexoptimset('cplex');
    o2 = options.cplex;
    n = fieldnames(o1);
    for i = 1:length(n)
        if isequal(o1.(n{i}),o2.(n{i}))
            options.cplex = rmfield(options.cplex,n{i});
        end
    end
catch
end