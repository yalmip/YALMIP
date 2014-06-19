function model = presolve_eliminatelinearratios(model)

% Try to linearize the model sum_j  aj xj/xi <= (==) b
% Common user mistake, can be written as sum_j  aj xj <= (==) b*xi
% (module signs etc)
% ex7_2_1 is a classical model exhibiting this poor modeling

if ~isempty(model.F_struc)
    for i = 1:model.K.l
        a = model.F_struc(i+model.K.f,2:end);
        b = model.F_struc(i+model.K.f,1);
        b0 = b;
        used = find(a);
        ok = 1;
        if all(model.variabletype(used)==4)
            % Signomial products
            ok = 1;
            k = [];
            for j = 1:length(used)
                m = model.monomtable(used(j),:);
                if nnz(m) <= 2
                    if nnz(m) == 1
                        % we have a term ai/xk
                        khere = find(m);
                        k = [k khere];
                        b = a(used(j));
                        a(used(j))=0;
                    else
                        % we have a term aixs/xk?
                        khere = find(m==-1);
                        k = [k khere];
                        s = find(m==1);
                        if isempty(s)
                            ok = 0;
                            break % nope, two signomials
                        end
                        a(s) = a(used(j));
                        a(used(j)) = 0;
                    end
                else
                    ok = 0;
                    break
                end
            end
            if ok
                if all(k==k(1))
                    if model.lb(k(1))>=0
                        a(k(1))=b0;
                        model.F_struc(i+model.K.f,:) = [b a];
                    elseif model.ub(k(1))<=0
                    end
                end
            end
        end
    end
end


