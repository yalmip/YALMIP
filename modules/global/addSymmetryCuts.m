function [p,p_lp] = addSymmetryCuts(p,p_lp)

for j = 1:length(p.sdpsymmetry)    
    if length(p.sdpsymmetry{j}.variables{1}) <= 4
        % We can enumerate easily infeasibles
        excludes = [];
        n = sqrt(size(p.sdpsymmetry{j}.dataBlock,1));
        combs = -dec2decbin(0:2^length(p.sdpsymmetry{j}.variables{1})-1,length(p.sdpsymmetry{j}.variables{1}))';
        for i = 1:size(combs,2)
            if min(eig(reshape(p.sdpsymmetry{j}.dataBlock*[1;combs(:,i)],n,n))) < -abs(p_lp.options.bnb.feastol)
                excludes = [excludes combs(:,i)];
            end
        end
        if ~isempty(excludes)
            newF = [];
            infeasible_combinations = unique(excludes','rows')';
            for k = 1:size(infeasible_combinations,2)
                % Local cut for reduced set
                [b,atemp] = exclusionCut(infeasible_combinations(:,k),-1);
                % Add that cut for every variable groups
                for s = 1:length(p.sdpsymmetry{j}.variables)
                    a = spalloc(1,length(p_lp.c),1);
                    a(p.sdpsymmetry{j}.variables{s}) = atemp;                    
                    newF = [newF;b a];                  
                end
            end
            p_lp = addInequality(p_lp,newF);
            % Delete, won't need these in the future
            p.sdpsymmetry{j}.dataBlock = [];
            p.sdpsymmetry{j}.variables = [];
        end
    end
end