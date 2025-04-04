function [f,q] = functionSeparation(p)

[monomtable,variabletype] = yalmip('monomtable');
evalvariables = yalmip('extvariables');
variabletype(evalvariables)=4;
q = p;
f = 0;
j = find(variabletype(p.lmi_variables(q.basis(1,2:end) ~= 0))>2);
f = f + sum(q.basis(1,1+j).*recover(p.lmi_variables(j(:)'))');
q.basis(1,1 + j) = 0;

j = find(variabletype(p.lmi_variables(q.basis(1,2:end) ~= 0))==2 | variabletype(p.lmi_variables(q.basis(1,2:end) ~= 0))==1);
for k = j(:)'
    s = find(monomtable(q.lmi_variables(k),:));
    if any(variabletype(s)>2)
        f = f + sum(q.basis(1,1+k).*recover(p.lmi_variables(k(:)'))');
        q.basis(1,1 + k) = 0;
    end
end


q = clean(q);
