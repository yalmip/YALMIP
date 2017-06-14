function Z = cleandoublefactors(Z)

if isempty(Z.midfactors)
    return
end

for i = 1:length(Z.midfactors)
    isnumeriq(i) = isnumeric(Z.midfactors{i});
end
isnumeriq = find(isnumeriq);
if length(isnumeriq)>1
    total = Z.leftfactors{isnumeriq(1)}*Z.midfactors{isnumeriq(1)}*Z.rightfactors{isnumeriq(1)};
    for i = isnumeriq(2:end)
        total = total + Z.leftfactors{i}*Z.midfactors{i}*Z.rightfactors{i};
    end
    allfactors = 1:length(Z.midfactors);
    keepfactors = setdiff(allfactors,isnumeriq(2:end));
    Z.midfactors{isnumeriq(1)} = total;
    Z.leftfactors{isnumeriq(1)} = speye(size(total,1));
    Z.rightfactors{isnumeriq(1)} = speye(size(total,2));
    Z.leftfactors = {Z.leftfactors{keepfactors}};
    Z.rightfactors = {Z.rightfactors{keepfactors}};
    Z.midfactors = {Z.midfactors{keepfactors}};
end
