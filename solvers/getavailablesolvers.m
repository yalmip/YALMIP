function [solvers,keep] = getavailablesolvers(findallsolvers,options);
    
solvers = definesolvers;
keep = ones(length(solvers),1);

% solverstocheck = 1:length(solvers);
% if length(options.solver)>0 & isempty(strfind(options.solver,'*')) & ~findallsolvers
%     for i = 1:length(solvers)
%         if strcmpi(options.solver,solvers(i).tag)
%            if ~solvers(i).usesother                              
%            end
%         end
%     end
% end

if ~findallsolvers
    for i = 1:length(solvers)
        isavailable = 1;
        j = 1;
        
        while (j <= length(solvers(i).checkfor)) & isavailable
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

if nargout == 1
    solvers = solvers(find(keep));
end
