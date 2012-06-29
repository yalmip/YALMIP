function model = mpt_de_mpmilp(Matrices,options,model)
                        
% Since we are recursively fixing binaries, big-M
% constraints will get active/inactive
[equalities,redundant] = mpt_detect_fixed_rows(Matrices);
if ~isempty(equalities)
    Matrices.Aeq = [Matrices.Aeq;Matrices.G(equalities,:)];
    Matrices.Beq = [Matrices.Beq;-Matrices.E(equalities,:)];
    Matrices.beq = [Matrices.beq;Matrices.W(equalities,:)];
    redundant = [redundant;equalities];
end
if ~isempty(redundant)
    Matrices.G(redundant,:) = [];
    Matrices.W(redundant,:) = [];
    Matrices.E(redundant,:) = [];
end

if all(Matrices.lb(Matrices.binary_var_index) == Matrices.ub(Matrices.binary_var_index))
    % Ok, we reached a node  
    model = mpt_solvenode(Matrices,Matrices.lb,Matrices.ub,Matrices,model,options);
else    
    % Find a non-fixed binary variable
    j = find(Matrices.lb(Matrices.binary_var_index) < Matrices.ub(Matrices.binary_var_index));
    j = Matrices.binary_var_index(j(1));

    % Solve up
    M = Matrices;
    M.lb(j) = 1;
    M.ub(j) = 1;
    if ~isempty(Matrices.SOSVariables)
        for i = 1:length(Matrices.SOS)
            if ismember(j, Matrices.binary_var_index(Matrices.SOS{i}))
                other = setdiff(Matrices.binary_var_index(Matrices.SOS{i}),j);
                M.lb(other) = 0;
                M.ub(other) = 0;
            end
        end
    end
                   
    % Simple bound tightening
    A = [Matrices.G -Matrices.E;Matrices.Aeq Matrices.Beq;-Matrices.Aeq -Matrices.Beq];
    b = [Matrices.W;Matrices.beq;Matrices.beq];
    [M.lb,M.ub,redundant,psstruct,infeasible] = tightenbounds(A,b,M.lb,M.ub,[],Matrices.binary_var_index,1:length(Matrices.lb));
           
    if any(M.lb(end-M.nx+1:end) >= M.ub(end-M.nx+1:end))      
        infeasible = 1;
    end
    if ~infeasible
        model = mpt_de_mpmilp(M,options,model);
    end

    % Solve down
    M = Matrices;
    M.lb(j) = 0;
    M.ub(j) = 0;
   %  [Matrices.lb(Matrices.binary_var_index) Matrices.ub(Matrices.binary_var_index)]
    [M.lb,M.ub,redundant,psstruct,infeasible] = tightenbounds(A,b,M.lb,M.ub,[],Matrices.binary_var_index,1:length(Matrices.lb));
    if any(M.lb(end-M.nx+1:end) >= M.ub(end-M.nx+1:end))         
        infeasible = 1;
    end
    if ~infeasible
        model = mpt_de_mpmilp(M,options,model);
    end
end
