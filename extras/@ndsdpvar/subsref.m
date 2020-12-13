function X = subsref(Y,refs)
% SUBSREF (overloaded)

for i = 1:length(refs.subs)
    if isa(refs.subs{i},'sdpvar')
        X = milpsubsref(Y,refs);
        return
    end
end
if isequal(refs.type,'()')
    if length(refs.subs)==3
        if length(Y.dim)==3 && isequal(Y.conicinfo,[-1 0])
            if length(refs.subs{1})==1 && length(refs.subs{2})==1 && ~isa(refs.subs{1},'char') && ~isa(refs.subs{2},'char') 
                if strcmp(refs.subs{3},':')
                    if refs.subs{1} >= 1 && refs.subs{1} <= Y.dim(1)
                        if refs.subs{2} >= 1 && refs.subs{2} <= Y.dim(2)
                            X = Y;
                            X.basis = [spalloc(X.dim(3),1,0) speye(X.dim(3))];
                            loc = 0:X.dim(1)*X.dim(2):X.dim(1)*X.dim(2)*(X.dim(3)-1);
                            i = refs.subs{1};
                            j = refs.subs{2};
                            X.lmi_variables = X.lmi_variables(loc + i + (j-1)*X.dim(1));
                            X.conicinfo = [0 0];
                            X.dim = [1 1 X.dim(3)];
                            return
                        end
                    end
                end
            end
        end
    elseif length(refs.subs)==1 && strcmp(refs.subs{1},':')
        X = Y;
        X.conicinfo = [0 0];
        X.dim = [prod(X.dim) 1];
        X = sdpvar([],[],[],[],[],[],[],struct(X));
        return
    end
end
X = Y;
base = reshape(1:size(Y.basis,1),X.dim);
base = subsref(base,refs);
X.basis = X.basis(base(:),:);
X.conicinfo = [0 0];
X.dim = size(base);
X = clean(X);