function p = presolve_sortrows(p)

if p.K.f > 0 
    % Sort the equalities such that equalities with few terms are at the
    % top. This will improve bound propagation from equalities. If a tie,
    % place equality with least number of variables with infinite bound on
    % top.
    s = p.lb-p.ub;
    b = p.F_struc(1:p.K.f,1);
    A = p.F_struc(1:p.K.f,2:end);
    [i,j] = sortrows([sum(A | A,2) sum(A(:,find(isinf(s))) | A(:,find(isinf(s))),2)],[1 2]);
    p.F_struc(1:p.K.f,:)=[];
    p.F_struc = [b(j) A(j,:);p.F_struc];
    if ~isempty(p.KCut.f)
        % Update the pointers to cut equalities
         [aux,newlocation] = ismember(p.KCut.f,j);
         p.KCut.f = newlocation;
    end
end