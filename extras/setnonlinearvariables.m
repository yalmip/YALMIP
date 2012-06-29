function x = setnonlinearvariables(p,x)

for i = 1:length(p.nonlinear)
    [row,pos,vals] = find(p.monomtable(p.nonlinear(i),:));
    x(p.nonlinear(i)) = prod((x(pos)').^(vals));
end

xevaled = x;%eros(1,length(p.c));
%xevaled(p.variabletype == 0) = x;
if ~isempty(p.evalVariables)
    % Experimental support for arbitrary functions
    if ~isempty(p.evalMap)
        for i = 1:length(p.evalMap)
            arguments = {p.evalMap{i}.fcn,xevaled(p.evalMap{i}.variableIndex)};
            arguments = {arguments{:},p.evalMap{i}.arg{2:end-1}};
            xevaled(p.evalVariables(i)) = feval(arguments{:});          
        end
    end
end

x = xevaled(:);
% for i = 1:length(p.nonlinear)
%     [row,pos,vals] = find(p.monomtable(p.nonlinear(i),:));
%     x(p.nonlinear(i)) = prod((x(pos)').^(vals));
% end
