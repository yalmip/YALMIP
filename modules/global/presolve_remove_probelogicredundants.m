function p = presolve_remove_probelogicredundants(p)
if nnz(p.Q) == 0 && isempty(p.evalMap) && p.K.l > 1 && ~isempty(p.binary_variables) && p.feasible
    probes = dec2decbin(0:7,3);
    ndeleted = 0;
    % How many rows is this variable involved in
    % We currently only check variables involved in 2
    % What we are looking for here are variables z and
    % models of the type
    %  x + z >= y
    %  1 + y >= x + z
    % This model has a z which makes it feasible for any (x,y),
    % so the constraints can effectively be removed
    % FIXME: Still have to recover z after solve!
    % FIXME: When constraints are removed, decrease variable count
    variablecount = sum(p.F_struc~=0,1);
    variablecount = variablecount(2:end);
    for i = p.binary_variables
        % Only in two constraints, and not in objective
        if variablecount(i)==2 && p.c(i)==0
            col = p.F_struc(:,1+i);
            [idx,~,~] = find(col);
            % and those are LP constraints
            if length(idx)==2 && all(idx > p.K.f) && all(idx <= p.K.f + p.K.l)
                bA = p.F_struc(idx,:);
                variables = find(any(bA(:,2:end),1));
                % There are in total three variables in these constraints
                if length(variables)==3 && all(p.isbinary(variables))
                    %  pos = sum(p.F_struc(:,1+variables) | p.F_struc(:,1+variables),1);
                    %  candidate = variables(pos==2);
                    % if ~isempty(candidate)
                    bA = bA(:,[1 1+variables]);
                    r = bA(:,1) + bA(:,2:end)*probes';
                    if all(any(r,1))
                        % This variable only occurs here
                        % and no matter what the other variables are
                        % it can be selected to make contraint feasible
                        % Hence it can be eliminated
                        % and constraints can be eliminated
                        p.F_struc(:,1+i) = 0;
                        p.F_struc(idx,[1 1+variables]) = 0;
                        ndeleted = ndeleted + 2;
                        % p.removableVariables = [p.removableVariables candidate];
                        p.removableVariables = [p.removableVariables i];
                    end
                    % end
                end
            end
        end
    end
    if p.options.verbose>1 && ndeleted > 0
        disp(['* Removed ' num2str(ndeleted) ' redundant logic constraints']);
        disp(['* Removed ' num2str(ndeleted/2) ' redundant binary variables']);
    end
end
p = presolve_empty_rows(p,0);