function [Fconv,no_changed,infeasible,Forigquad] = convertquadratics(F)
%CONVERTQUADRATICS Internal function to extract quadratic constraints

% ******************************
% LINEAR?
% ******************************

infeasible = 0;
Fconv = F;
no_changed = 0;
Forigquad = [];

if islinear(F)
    return
end

[monomtable,variabletype] = yalmip('monomtable');

if any(variabletype == 4)
    if issigmonial(F)
        return
    end
end

Fconv = lmi;
no_changed = 0;
i_changed = [];
% Make sure bounds a pre-processed, in case we look for rotated cone which
% requires us to know lower bounds
nv = yalmip('nvars');
LU = yalmip('getbounds',1:nv);
LUhere = getbounds(F,[],LU);
for i = 1:1:length(F)
    if max(variabletype(getvariables(F(i)))) <= 1
        % Definitely no quadratic to model as all variables are bilinear at
        % most
         Fconv = Fconv + F(i);
    elseif is(F(i),'element-wise') & ~is(F(i),'linear') & ~is(F(i),'sigmonial')
        % f-c'*x-x'*Q*x>0
        fi = sdpvar(F(i));fi = fi(:);
        %[Qs,cs,fs,x,info] = vecquaddecomp(fi);
        for j = 1:length(fi)
            fij = fi(j);
            if isa(fij,'double')
                if fij < 0
                    infeasible = 1;                    
                    warning('The problem is trivially infeasible. You have a term POSITIVE NUMBER < 0.');
                    return
                end
            end
           % Q = Qs{j};c = cs{j};f = fs{j};
           if max(variabletype(getvariables(fij))) <= 1
               % There are at most bilinear terms, so it cannot be convex
               info = 1;
           else
               [Q,c,f,x,info] = quaddecomp(fij);
           end
            if info==0
                if nnz(Q)==0
                    % Oh, linear,...
                    Fconv = Fconv + (fi(j)>=0);
                else
                    % Yes, quadratic, but convex?
                    % Change sign definitions
                    Q = -Q;
                    c = -c;
                    f = -f;
                    % Semi-definite case when only part of x in Q
                    % Occurs, e.g, in constraints like y'*Q*y < t
                    used = find(any(Q));Qred=Q(:,used);Qred = Qred(used,:);xred = x(used);
                    [R,p]=chol(Qred);
                    done = 0;
                    if p
                        % Safety check to account for low rank problems
                        if all(eig(full(Qred))>=-1e-12)
                            [u,s,v]=svd(full(Qred));
                            r=find(diag(s)>1e-12);
                            R=(u(:,r)*sqrt(s(r,r)))';
                            R(abs(R)<eps)=0;
                            p=0;
                        else
                            % Try to detect rotated SOCP (x-d)'*A*(x-d)+k<= C*y*z
                            yzCandidates = find(~diag(Qred));
                            if length(yzCandidates) == 2 && nnz(c(yzCandidates))==0
                                %LU = yalmip('getbounds',getvariables(xred(yzCandidates)));
                                LU = LUhere(getvariables(xred(yzCandidates)),:);
                                if all(LU(:,1)>=0)
                                    yzSubQ = Qred(yzCandidates,yzCandidates);
                                    if yzSubQ(1,2) < 0
                                        C = -2*yzSubQ(1,2);
                                        xCandidates = setdiff(1:length(used),yzCandidates);
                                        A = Qred(xCandidates,xCandidates);
                                        [B,p]=chol(A);
                                        if ~p
                                            y = xred(yzCandidates(1));
                                            z = xred(yzCandidates(2));
                                            d = -A\(c(xCandidates)/2);
                                            k = f-d'*A*d;
                                            if abs(k)<= 1e-12
                                                Fconv=Fconv + lmi(rcone(B*(xred(xCandidates)-d),.5*C*y,z));
                                                no_changed = no_changed + 1;
                                                i_changed = [i_changed i];
                                                done = 1;
                                            elseif k >= -1e-12
                                                Fconv=Fconv + lmi(rcone([B*(xred(xCandidates)-d);sqrt(abs(k))],.5*C*y,z));
                                                no_changed = no_changed + 1;
                                                i_changed = [i_changed i];
                                                done = 1;
                                            else
                                                p = 1;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if p==0 && ~done
                        % Write as second order cone
                        d = -c'*x-f;
                        if isa(d,'double')
                            if d<0
                                infeasible = 1;
                                return
                            else
                                Fconv=Fconv + lmi(cone([R*xred],sqrt(d)));
                            end
                        else
                            if length(c) == length(xred)
                                ctilde = -(R')\(c/2);                             
                                Fconv=Fconv + lmi(cone([R*xred;.5*(1-d)],.5*(1+d)));                             
                            else
                                Fconv=Fconv + lmi(cone([R*xred;.5*(1-d)],.5*(1+d)));
                            end
                        end
                        no_changed = no_changed + 1;
                        i_changed = [i_changed i];                        
                    elseif ~done
                        Fconv = Fconv + lmi(fi(j));
                    end
                end
            else
                Fconv = Fconv + lmi(fi(j));
            end
        end
    else
        Fconv = Fconv + F(i);
    end
end
if ~isempty(i_changed)
    Forigquad = F(i_changed);
else
    Fconv = F;
end
