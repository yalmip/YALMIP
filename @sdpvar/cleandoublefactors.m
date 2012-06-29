function Z = cleandoublefactors(Z)

if isempty(Z.midfactors)
    return
end

for i = 1:length(Z.midfactors)
    isdouble(i) = isa(Z.midfactors{i},'double');
end
isdouble = find(isdouble);
if length(isdouble)>1
    total = Z.leftfactors{isdouble(1)}*Z.midfactors{isdouble(1)}*Z.rightfactors{isdouble(1)};
    for i = isdouble(2:end)
        total = total + Z.leftfactors{i}*Z.midfactors{i}*Z.rightfactors{i};
    end
    allfactors = 1:length(Z.midfactors);
    keepfactors = setdiff(allfactors,isdouble(2:end));
    Z.midfactors{isdouble(1)} = total;
    Z.leftfactors{isdouble(1)} = eye(size(total,1));
    Z.rightfactors{isdouble(1)} = eye(size(total,2));
    Z.leftfactors = {Z.leftfactors{keepfactors}};
    Z.rightfactors = {Z.rightfactors{keepfactors}};
    Z.midfactors = {Z.midfactors{keepfactors}};
end
