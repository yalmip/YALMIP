function mu = meanFromDistribution(distribution)


if isa(distribution.name, 'function_handle') && strcmpi(func2str(distribution.name),'random')
    switch  distribution.parameters{1}
        case 'normal'
            mu = distribution.parameters{2};
            
        case 'exponential'
        otherwise
    end
end

