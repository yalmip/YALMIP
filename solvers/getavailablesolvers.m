function [solvers,keep,allsolvers] = getavailablesolvers(findallsolvers,options);
    
solvers = definesolvers;allsolvers = solvers;
keep = ones(length(solvers),1);

if ~findallsolvers
    for i = 1:length(solvers)
        isavailable = 1;
        j = 1;
        
        while (j <= length(solvers(i).checkfor)) && isavailable
            s = exist(solvers(i).checkfor{j},'file');
            s = (s~=0) & (s~=7);
            isavailable = isavailable & s;
            j = j + 1;
        end
        if ~isavailable
            keep(i)=0;
        end
    end
end

if nargout == 1 || nargout==3
    solvers = solvers(find(keep));
end
