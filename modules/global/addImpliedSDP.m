function p = addImpliedSDP(p)
if p.K.q(1)==0 && p.K.s(1) > 0
    % Search for rows where there is a constant, and a diagonal term 
    % which thus has to be structly positive, we might leads us to a
    % constraints of the type sum x_i >= 1, if only non-negative integer
    % variables are involved, entering with positive coefficient (and v.s)
    newF = [];
    top = sum(p.K.q) + p.K.f + p.K.l + 1;
    for i = 1:length(p.K.s)
        F = p.F_struc(top:top + p.K.s(i)^2-1,:);
        % Constant term in SDP
        X0 = reshape(F(:,1),p.K.s(i),p.K.s(i));
        % linear terms 
        XX = reshape(sum(F(:,2:end) ~= 0,2),p.K.s(i),p.K.s(i));
        candidates = find(any(X0-diag(diag(X0)),2) & (diag(X0)<=0) & ~any(XX-diag(diag(XX)),2));
        if ~isempty(candidates);
            I = speye(p.K.s(i),p.K.s(i));
            pos = find(I(:));
            for j = candidates(1:end)'
                row = F(pos(j),2:end);
                used = find(row);
                if ~isempty(used)
                    if ~any(row(p.noninteger_variables))
                        if all(row <=0) && all(p.ub(used)<=0)
                            row(find(row)) = 1;
                            newF = [newF;-1 -row];
                        elseif all(row >=0) && all(p.lb(used)>=0)
                            row(find(row)) = 1;
                            newF = [newF;-1 row];
                        end
                    end
                end
            end
        end
        top = top + p.K.s(i)^2;
    end
    if size(newF,1)>0
        p.F_struc = [p.F_struc(1:p.K.f + p.K.l,:);
            newF;
            p.F_struc(1+p.K.f+p.K.l:end,:)];
        p.K.l = p.K.l + size(newF,1);
    end
    
    % Search for trivial variables entering in diagonal elements alone, and
    % not entering anywhere else in model. Bad modelling?
    % First, they should not appear in LP cone (can be generalized by
    % checking that it acts in the correct direction)
    % TODO: Check for possibility of different signs etc, currently only
    % implemented for typical test cases
    candidates = find(~any(p.F_struc(1:p.K.f + p.K.l,2:end),1));
    fixable = nan(1,length(candidates));   
    
    top = sum(p.K.q) + p.K.f + p.K.l + 1;
    for j = 1:length(p.K.s)
        F = p.F_struc(top:top + p.K.s(j)^2-1,:);
        for i = 1:length(candidates(:)')
            X0 = reshape(F(:,1 + candidates(i)),p.K.s(j),p.K.s(j));
            % Diagonal matrix?
            d = diag(X0);
            if nnz((X0-diag(d)))==0                
                if all(sign(d)<=0) && p.c(i)>=0 && (isnan(fixable(i)) || fixable(i)==-1)
                    fixable(i) = -1;
                elseif all(sign(d)>=0) && p.c(i)<=0 && (isnan(fixable(i)) || fixable(i)==1)
                    fixable(i) = 1;
                else
                    fixable(i) = 0;
                end
            else
                fixable(i) = 0;
            end          
        end
        top = top + p.K.s(j)^2;
    end
    p.adjustable = candidates(find(fixable));
    for i = candidates(find(~isnan(fixable)))
%         if fixable(i) == -1
%             if ~isinf(p.lb(i))
%                 %             p.ub(i) = p.lb(i);
%                 %             p.F_struc(:,1) = p.F_struc(:,1) + p.F_struc(:,i+1)*p.lb(i);
%                 %             p.F_struc(:,i+1) = 0;
%                 %         else
%                 %             % TODO whole row/column should be removed from the LMI!
%             end
%         elseif fixable(i) == 1
%             if ~isinf(p.ub(i))
%                 %             p.ub(i) = p.lb(i);
%                 %             p.F_struc(:,1) = p.F_struc(:,1) + p.F_struc(:,i+1)*p.lb(i);
%                 %             p.F_struc(:,i+1) = 0;
%                 %         else
%                 %             % TODO whole row/column should be removed from the LMI!
%             end
%         end
        
    end
end

