function LU = extract_bounds_from_milpsubsref_operator(LU,extstruct,extvariables,i)
arg = extstruct(i).arg{1};
epi = getvariables(extstruct(i).var);
[M,m] = derivebounds(reshape(arg,[],1),LU);
LU(epi,1) = min(m);
LU(epi,2) = max(M);
for j = 1:length(extstruct(i).arg{2}.subs)
    index = extstruct(i).arg{2}.subs{j};
    if isa(index,'sdpvar')
        if isequal(getbase(index),[0 1])
            LU(getvariables(index),1) = max(LU(getvariables(index),1),1);
            LU(getvariables(index),2) = min(LU(getvariables(index),2),numel(arg));
        end
    end
end
