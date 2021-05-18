function [p,lower,upper] = find_simple_variable_bounds(p);
% This internal function extracts variable bounds for variables that only
% have upper and lower bounds, and are not involved in any other
% constraints. If the variable is involved in any other constraint, the
% bound is set to inf

if any(p.K.q) | any(p.K.s)
    lower = -inf(length(p.c),1);
    upper = inf(length(p.c),1);
else
    [lower,upper,used_rows_eq,used_rows_lp] = find_lp_bounds(p.F_struc,p.K);
    used_rows_eq = used_rows_eq(~any(full(p.F_struc(used_rows_eq,1+find(p.variabletype~=0))),2)); 
    used_rows_lp = used_rows_lp(~any(full(p.F_struc(p.K.f + used_rows_lp,1+find(p.variabletype~=0))),2));
    if ~isempty(used_rows_lp)
        p_temp = p;
        p_temp.F_struc(p.K.f+used_rows_lp,:)=[];
        p_temp.K.l = p.K.l - length(used_rows_lp);

        if size(p.F_struc,1) > 0
            % These variables are still used in some other constraints
            still_used = find(sum(abs(p_temp.F_struc(:,2:end)),1) > 0);
            if ~isempty(still_used)
                % we have to keep these variables
                lower(still_used) = -inf;
                upper(still_used) = inf;            
                
                keep_rows_eq = find(sum(abs(p.F_struc(1:p.K.f,1+still_used)),2) > 0);
                keep_rows_lp = find(sum(abs(p.F_struc(1+p.K.f:p.K.f+p.K.l,1+still_used)),2) > 0);
                keep_rows_c = find(sum(abs(p.F_struc(1+p.K.f+p.K.l:end,1+still_used)),2) > 0);
%                 if ~isempty(keep_rows_c)% any(keep_rows>(p.K.l + p.K.f))
%                     error('Tell johan to fix the SDP case in find_simple_variable_bounds!')
%                 end
                used_rows_eq = setdiff(used_rows_eq,keep_rows_eq);
                used_rows_lp = setdiff(used_rows_lp,keep_rows_lp);
                p.F_struc(used_rows_lp+p.K.f,:)=[];
                p.F_struc(used_rows_eq,:)=[];
                p.K.l = p.K.l - length(used_rows_lp);                
                p.K.f = p.K.f - length(used_rows_eq);
                                  
%                 keep_rows = find(sum(abs(p.F_struc(:,1+still_used)),2) > 0);
%                 if any(keep_rows>(p.K.l + p.K.f))
%                       error('Tell johan to fix the SDP case in find_simple_variable_bounds!')
%                 end
%                 used_rows = used_rows + p.K.f;
%                 used_rows = setdiff(used_rows,keep_rows);
%                 p.F_struc(used_rows,:)=[];
%                 p.K.l = p.K.l - nnz(used_rows>p.K.f);
%                 p.K.f = p.K.f - nnz(used_rows<=p.K.f);
            else
                p = p_temp;
            end
        else
            p = p_temp;
        end
    end
end

