function x = evaluate_nonlinear(p,x,qq)

% FIX: We have to apply computations to make sure we are evaluating
% expressions such as log(1+sin(x.^2).^2) correctly

if ~isempty(p.bilinears) & all(p.variabletype <= 2) & length(p.evalMap)==0
    x(p.bilinears(:,1)) = x(p.bilinears(:,2)).*x(p.bilinears(:,3));
else
    oldx = 0*p.c;old(1:length(x))=x;
    %try
    x = process_polynomial(x,p);
    %catch
    %    1
    %end
    x = process_evaluations(x,p);
    while norm(x - oldx)>1e-8
        oldx = x;
        x = process_polynomial(x,p);
        x = process_evaluations(x,p);
    end
end

function x = process_evaluations(x,p)
for i = 1:length(p.evalMap)
    arguments = {p.evalMap{i}.fcn,x(p.evalMap{i}.variableIndex)};
    arguments = {arguments{:},p.evalMap{i}.arg{2:end-1}};
    x(p.evalVariables(i)) = feval(arguments{:});
    if ~isempty(p.bilinears)
        x = process_bilinear(x,p);
    end
end

function x = process_bilinear(x,p)
x(p.bilinears(:,1)) = x(p.bilinears(:,2)).*x(p.bilinears(:,3));

function x = process_polynomial(x,p)
x = x(1:length(p.c));
nonlinear = find(p.variabletype);
x(nonlinear) = prod(repmat(x(:)',length(nonlinear),1).^p.monomtable(nonlinear,:),2);
