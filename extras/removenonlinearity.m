function  p = removenonlinearity(p)
p.variabletype = 0*p.variabletype;
p.evalMap = [];
if ~isempty(p.x0)
    if nnz(p.x0) == 0
        p.x0 = zeros(length(p.c),1);
    end
else
    % FIXME: Initials should be made consistent
    p.x0 = zeros(length(p.c),1);
end