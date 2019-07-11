function p = compile_bilinearslist(p);
% Precomputed list of bilinear expressions, used in when evaluating value
% which is done in heuristics code
Bilinears = find(p.variabletype==1);
BilinearsList = zeros(length(p.c),2);
for i = Bilinears
    vars = find(p.monomtable(i,:));
    BilinearsList(i,:) = vars(:)';
end
p.Bilinears = Bilinears;
p.BilinearsList = BilinearsList;