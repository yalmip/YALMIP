function [p,p_lp] = addSymmetryCuts(p,p_lp)

for j = 1:length(p.sdpsymmetry)    
    if length(p.sdpsymmetry{j}.variables{1}) <= 4
        % We can enumerate easily infeasibles
        excludes = [];
        n = sqrt(size(p.sdpsymmetry{j}.dataBlock,1));
        vars = p.sdpsymmetry{j}.variables{1};
        % Now enumerate, for historical reasons
        % also accept {-1,0}-models so translate and scale
        L = p.lb(vars);
        U = p.ub(vars);
        m = length(vars);
        combs = dec2decbin(0:2^m-1,m)';        
        combs = repmat(L,1,size(combs,2)) + repmat(U-L,1,size(combs,2)).*combs;
        combs = unique(combs','rows')';
        feasible = ones(1,size(combs,2));
        for i = 1:size(combs,2)
            if min(eig(reshape(p.sdpsymmetry{j}.dataBlock*[1;combs(:,i)],n,n))) < -abs(p_lp.options.bnb.feastol)
                excludes = [excludes combs(:,i)];
                feasible(i) = 0;
            end
        end
        % We can model the inconsistencies in various ways       
        % 1. Binary variables sum less
        % 2. Binary variables sum larger
        % 3. Generic, just exclude                           
        if ~isempty(excludes)
            done = 0;
            sum_feas = sum(combs(:,find(feasible)),1);
            sum_infeas = sum(excludes,1); 
            newF = [];
            for k = min(sum(combs,1)):max(sum(combs,1))
                if all(sum_feas<=k) && all(sum_infeas>k)
                    for s = 1:length(p.sdpsymmetry{j}.variables)
                        a = spalloc(1,length(p_lp.c),1);
                        a(p.sdpsymmetry{j}.variables{s}) = -1;
                        newF = [newF;k a];                          
                    end                    
                    done = 1;
                    break
                elseif  all(sum_feas>=k) && all(sum_feas<k)
                    for s = 1:length(p.sdpsymmetry{j}.variables)
                        a = spalloc(1,length(p_lp.c),1);
                        a(p.sdpsymmetry{j}.variables{s}) = 1;
                        newF = [newF;-k a];                       
                    end 
                    done = 1;
                    break
                end
            end
            if ~done   
                % No good structure, just add no-goods
                infeasible_combinations = unique(excludes','rows')';
                for k = 1:size(infeasible_combinations,2)
                    % Local cut for reduced set
                    if min(min(excludes))==-1
                        binary_type = -1;
                    else
                        binary_type = 1;
                    end
                    [b,atemp] = exclusionCut(infeasible_combinations(:,k),binary_type);
                    % Add that cut for every variable groups
                    for s = 1:length(p.sdpsymmetry{j}.variables)
                        a = spalloc(1,length(p_lp.c),1);
                        a(p.sdpsymmetry{j}.variables{s}) = atemp;
                        newF = [newF;b a];
                    end
                end
            end
            p_lp = addInequality(p_lp,newF);
            % Delete, won't need these in the future
            p.sdpsymmetry{j}.dataBlock = [];
            p.sdpsymmetry{j}.variables = [];
        end
    end
end