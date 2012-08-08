function featureSupported = solverCapable(solverList,solverName,feature);

i = 1;
featureSupported = 1;
while i <= length(solverList)
    if strcmpi(solverList(i).tag,solverName)
        % dynamic fields not supported in old versions...
        featureSupported = eval(['solverList(i).' feature]);
        return
    end
    i = i+1;
end
            