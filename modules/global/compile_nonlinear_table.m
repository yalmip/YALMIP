function p = compile_nonlinear_table(p)
linears = find(p.variabletype == 0);
nonlinears =  find(p.variabletype > 0);
bilinears   = [];
for i = 1:length(nonlinears)
    if p.variabletype(nonlinears(i))<3
        z = find(p.monomtable(nonlinears(i),:));
        if length(z)==1
            bilinears = [bilinears;nonlinears(i) z z];
        else
            bilinears = [bilinears;nonlinears(i) z(1) z(2)];
        end
    end
end
nonlinears = union(nonlinears,p.evalVariables);
linears = setdiff(linears,p.evalVariables);

p.linears = linears;
p.bilinears = bilinears;
p.nonlinears = nonlinears;

Quadratics = find(p.variabletype==2);
QuadraticsList = zeros(length(p.c),2);
for i = Quadratics
    vars = find(p.monomtable(i,:));
    QuadraticsList(i,:) = vars(:)';
end
p.Quadratics = Quadratics;
p.QuadraticsList = QuadraticsList;
