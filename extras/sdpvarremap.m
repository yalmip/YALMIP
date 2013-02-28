function X = sdpvarremap(X,old,new)

[index,pos] = ismember(old, getvariables(X));
index = find(index);
pos = pos(index);
if ~isempty(pos)
    temp = recover(getvariables(X));
    temp(pos) = new(index);
    temp = getbase(X)*[1;temp];
    X = reshape(temp,size(X));
end