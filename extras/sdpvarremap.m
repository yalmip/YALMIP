function X = sdpvarremap(X,old,new)

[index,pos] = ismember(old, getvariables(X));
if nnz(index)>0
    index = find(index);
    pos = pos(index);
   
     B = getbase(X);
    if 1
       temp = getvariables(X);
           
       cnew = temp(pos);
       nopos = setdiff(1:length(temp),pos);
       
       if isempty(nopos)
           temp = B(:,1) + B(:,1+pos)*new(index);
       else
           temp = B(:,1) + B(:,1+pos)*new(index)+B(:,1+nopos)*recover(temp(nopos));
       end
    elseif 0
        variables = getvariables(X);
        keeppos = setdiff(variables,index)        
               
       
        if ~isempty(keeppos)
            keepvariables = variables(keeppos);
            temp = B(:,1)+B(:,keeppos + 1)*recover(keepvariables)+B(:,pos + 1)*new(index);
        else
            temp = B(:,1)+B(:,pos + 1)*new(index);
        end
    else
        temp = recover(getvariables(X));
        temp(pos) = new(index);
        temp = getbase(X)*[1;temp];
    end
    X = reshape(temp,size(X));
end