function p = presolve_downforce(p)

p.downForce = {};
if p.K.l == 0
    return
end
% Extract linear inequalities
ff = p.F_struc(startofLPCone(p.K):startofLPCone(p.K)+p.K.l-1,:);
% Search for y >= x1 + x2 +...
candidates = find(ff(:,1)==0);
for i = 1:length(candidates)
    row = p.F_struc(candidates(i),2:end);
    xx = find(row==-1);
    yy = find(row==1);    
    if length(yy)==1 && length(xx)+1==nnz(row)
        % Yes, this matches y>=x1+x2+x3+...
        if all(ismember(find(row),p.binary_variables))
            p.downForce{end+1}.forcing = yy;
            p.downForce{end}.forced = xx;            
        end
    end
end
