function p = reduce_bilinear_branching_variables(p)
if ~p.originallyNonlinearConstraints && p.solver.lowersolver.objective.quadratic.convex & p.problemclass.objective.quadratic.nonconvex==0 && isempty(p.evalMap)
    % Setup quadratic
    Q_ = p.Q;
    for i = 1:size(p.bilinears,1)
        if p.c(p.bilinears(i,1))
            Q_(p.bilinears(i,2),p.bilinears(i,3)) = p.c(p.bilinears(i,1))/2;
            Q_(p.bilinears(i,2),p.bilinears(i,3)) = Q_(p.bilinears(i,3),p.bilinears(i,2))+p.c(p.bilinears(i,1))/2;
        end
    end
    if nnz(Q_)>0
        % Cheap & quick vs expensive and correct
        if ((nnz(diag(Q_))==nnz(Q_)) & all(diag(Q_)>=0)) | all(eig(full(Q_))>-1e-12)
            Used_in_F = find(any(p.F_struc(:,2:end),1));
            Used_in_F = intersect(Used_in_F,p.bilinears(:,1));
            %  old_branch =  p.branch_variables;
            p.branch_variables = [];
            for i = 1:size(p.bilinears,1)
                j = p.bilinears(i,1);
                if ismember(j,Used_in_F)
                    p.branch_variables = [p.branch_variables p.bilinears(i,2:3)];
                end
            end
            p.branch_variables = unique(p.branch_variables);
            %   p.branch_variables = intersect(old_branch,p.branch_variables);
        end
    end
end

% Do not branch in auxilliary variables introduced to simply normalize the
% nonlinear operators, such as cos(2x + y) replaced with cos(z), z==2x+y
% if ~isempty(p.evalMap)
%     if any(p.K.f)
%         allInArg = [];
%         for i = 1:length(p.evalMap)
%             allInArg = [allInArg p.evalMap{i}.variableIndex];
%         end
%         rmv = [];
%         for i = 1:length(p.evalMap)
%             inArg = p.evalMap{i}.variableIndex;
%             if length(inArg) == 1
%                 coeffs = p.F_struc(:,1 + inArg);
%                 pos = find(coeffs);
%                 if length(pos) == 1 & pos<=p.K.f
%                     watch out for x1+x2==1, minimize sin(x1)+sin(x2)
%                     a = p.F_struc(pos,2:end);
%                     a(inArg)=0;
%                     if ~any(ismember(find(a),allInArg))
%                         rmv = [rmv inArg];
%                     end
%                 end
%             end
%         end
%         p.branch_variables = setdiff( p.branch_variables,rmv);
%     end
% end