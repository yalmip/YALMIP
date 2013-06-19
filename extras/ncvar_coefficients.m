function [c,v] = ncvar_coefficients(p,x)

nonCommutingTable         = yalmip('nonCommutingTable');
[monomtable,variabletype] = yalmip('monomtable');
if size(monomtable,1)>size(nonCommutingTable,1)
    nonCommutingTable((1+size(nonCommutingTable,1)):(size(monomtable,1)),1) = (1+size(nonCommutingTable,1)):(size(monomtable,1));
end

base = getbase(p);
vars = getvariables(p);
xvar = getvariables(x);
c = [];
v = [];
vvars = [];
for k = 1:length(vars)
    i = vars(k);
    monTerms = nonCommutingTable(i,:);
    if isnan(monTerms(1))
        % pure non-commuting, hence it can not involve x
        c = [c; recover(i)];
        v = [v; 1];
        vvars = [vvars;0];
    else
        % [v(x) non-commuting]. Call recursively on v(x)
        [ctemp,vtemp] = coefficients(recover(monTerms(1)),x);
        if ~isempty(ctemp)
            multiplier = 1;
            for j = 2:max(find(monTerms))
                h = recover(monTerms(j));
                multiplier = multiplier*h;
            end
            c = [c;ctemp*multiplier];
            v = [v;vtemp];
            if isa(vtemp,'double')
                vvars = [vvars;0];
            else                
                vvars = [vvars;getvariables(vtemp)];
            end
        end
    end
end

% Compress to unique terms in basis
% vvar = [];
% for i = 1:length(v)
%     vvar = [vvar getvariables(v(i))];
% end
[unique_monomials,keep,j] = unique(vvars);
M = spalloc(length(unique_monomials),length(v),0);
for i = 1:length(unique_monomials)
    M(i,unique_monomials(i) == vvars) = 1;
end
c = M*c;
v = v(keep);