function model = compile_quadraticslist(model)
Quadratics = find(model.variabletype==2);
QuadraticsList = zeros(length(model.c),2);
for i = Quadratics
    vars = find(model.monomtable(i,:));
    QuadraticsList(i,:) = vars(:)';
end
model.Quadratics = Quadratics;
model.QuadraticsList = QuadraticsList;