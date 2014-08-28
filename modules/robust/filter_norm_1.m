function [F,feasible] = filter_norm_1(all_f,all_c_w,all_c_x,all_Q_xw,x,Zmodel,allw,X,Q_xx,VariableType)

feasible = 1;
% As a first step, we figure out the radius
r = Zmodel.r;
center = Zmodel.center;
F = ([]);
ci_basis = all_c_w';
lastBici = [];
for i = 1:length(all_f)
    Bi = 2*all_Q_xw(:,length(x)*(i-1)+1:length(x)*i)';
    bi = all_c_x(length(x)*(i-1)+1:length(x)*i);
    if (nnz(ci_basis(:,i))==0) & nnz(Bi)==0
        F = F + (X(i)>=0);        
    else
        ci = ci_basis(:,i);
        di = all_f(i);
        
        % Shift |w-center|, wtilde = w-center i.e. w=wtilde+center
        di = di + ci(1:length(center),:)'*center;
        bi = bi + Bi(:,1:length(center))*center;
        
        used = find(full(any([full(Bi') ci],2)));
        %ci = ci(used);
        %Bi = Bi(:,used);
        % Shift |w-center|, wtilde = w-center i.e. w=wtilde+center
        %di = di + ci'*center;
        %bi = bi + Bi*center;
        if isequal(lastBici,[Bi' ci])
            F = F + (x'*Q_xx{i}*x+bi'*x + di - r*s >= 0);% + (-s<Bi'*x+ci<s);
        else
            if nnz(Bi)==0
                s = norm(full(ci),inf);
            else
                s = sdpvar(1,1);
                F =  F + (-s<=Bi'*x+ci<=s);
            end
            F = F + (x'*Q_xx{i}*x+bi'*x + di - r*s >= 0);
            lastBici = [Bi' ci];
        end
    end
end