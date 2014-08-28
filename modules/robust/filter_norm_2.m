function [F,feasible] = filter_norm_2(all_f,all_c_w,all_c_x,all_Q_xw,x,Zmodel,allw,X,Q_xx,VariableType)

feasible = 1;
% As a first step, we figure out the radius
r = Zmodel.r;
center = Zmodel.center;
F = ([]);
ci_basis = all_c_w';
lastBici = [];
lastusedrows = [];
% (bi' + (Bi*w)')*x + (ci'*w + di).
for i = 1:length(all_f)
    Bi = 2*all_Q_xw(:,length(x)*(i-1)+1:length(x)*i)';
    bi = all_c_x(length(x)*(i-1)+1:length(x)*i);
    if (nnz(ci_basis(:,i))==0) & nnz(Bi)==0
        F = F + (X(i)>=0);
        % This constraint row is constant
    else
        ci = ci_basis(:,i);
        di = all_f(i);
        used = find(full(any([full(Bi') ci],2)));
        %        ci = ci(used);
        %        Bi = Bi(:,used);
        % Shift |w-center|, wtilde = w-center i.e. w=wtilde+center
        di = di + ci'*center;
        bi = bi + Bi*center;
        if isequal(abs(lastBici),abs([Bi' ci]))            
             F = F + (x'*Q_xx{i}*x+bi'*x + di - r*s >= 0);            
      %       F = F + (x'*Q_xx{i}*x+bi'*x + di - r*norm(full(Bi(used_rows,used)))*s >= 0);            
        else
            s = sdpvar(1,1);
            if length(used)==1
                if nnz(Bi(:,used))==0
                    s = norm(full(ci(used)));
                    F = F + (x'*Q_xx{i}*x+bi'*x + di - r*s >= 0);
                else
                    F = F + (x'*Q_xx{i}*x+bi'*x + di - r*s >= 0) + (-s<=Bi(:,used)'*x+ci(used)<=s);                   
                end
                used_rows = 1:size(Bi,1);
                lastBici = [Bi' ci];
            else
                used_rows = find(any(full(Bi(:,used)')));
                if length(used_rows)==1 & nnz(ci(used))==0
                    % Special case norm(double vector*scalar)
                    if isequal(used_rows,lastusedrows)
                        F = F + (x'*Q_xx{i}*x+bi'*x + di - r*norm(full(Bi(used_rows,used)))*lasts >= 0)  ;
                    else
                        F = F + (x'*Q_xx{i}*x+bi'*x + di - r*norm(full(Bi(used_rows,used)))*s >= 0) + (-s <= x(used_rows)<=s);
                        lastBici = [Bi' ci];
                        lastusedrows = used_rows;
                        lasts = s;
                    end
                else
                    if nnz(Bi)==0
                        s = norm(full(ci(used)));
                        F = F + (x'*Q_xx{i}*x+bi'*x + di - r*s >= 0) ;
                    else
                        F = F + (x'*Q_xx{i}*x+bi'*x + di - r*s >= 0) + (cone(Bi(:,used)'*x+ci(used),s));
                    end
                    lastBici = [Bi' ci];
                end
            end
        end
    end
end