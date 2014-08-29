function Y=simplify(X,how)
%SIMPLIFY  Reduce PWA complexity

if nargin <2
    how = 'greedy';
end

variables = getvariables(X);
extstruct = yalmip('extstruct',variables(1));
if ~isempty(extstruct)
    if isequal(extstruct.fcn,'pwa_yalmip')
        extstruct.arg{1}{1}.Fi = extstruct.arg{1}{1}.Bi;
        extstruct.arg{1}{1}.Gi = extstruct.arg{1}{1}.Ci;
        simplified = mpt_simplify(extstruct.arg{1}{1},how);
        simplified.Bi = simplified.Fi;
        simplified.Ci = simplified.Gi;
        Y = pwf(simplified,extstruct.arg{2},extstruct.arg{3});
    end
end
