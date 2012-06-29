function p = eliminateBinary(p,binaries)

vars = p.lmi_variables;
[mt,vt] = yalmip('monomtable');


mt = mt(vars,:);
mt(:,binaries) = min(mt(:,binaries),1);%rem(mt(:,binaries),2);

used_variables = find(any(mt,1));
x = recover(used_variables)';
new_monoms = [];
mt = mt(:,used_variables);

y = recovermonoms(mt,x);
y = p.basis*[1;y];
if isa(y,'sdpvar')
    % copy data to p
    p.basis = y.basis;
    p.lmi_variables = y.lmi_variables;
else
    p = y;
end
