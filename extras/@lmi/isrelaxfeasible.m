function res = isrelaxfeasible(F)


% Author Johan Löfberg
% $Id: isrelaxfeasible.m,v 1.3 2005-02-04 10:10:27 johanl Exp $

% Check if solution avaliable
currsol = evalin('caller','sdpvar(''getSolution'')');
if isempty(currsol)
    disp('No solution available.')
    return
end

nlmi = size(F.clauses,2);
spaces = ['                                    '];
if (nlmi == 0) & (neq == 0)
    disp('empty LMI')
    return
end

lmiinfo{1} = 'LMI';
lmiinfo{2} = 'Element-wise';
lmiinfo{3} = 'Equality constraint';
lmiinfo{4} = 'Second order cone constraint';
lmiinfo{5} = 'Rotated Lorentz constraint';

header = {'ID','Constraint','Type','Residual (should be > 0)','Tag'};

if nlmi>0
    for j = 1:nlmi
        F0 = relaxdouble(F.clauses{j}.data);
        if any(isnan(F0(:)))
            res = NaN;
        else
            switch F.clauses{j}.type
                case 1
                    res = min(eig(F0));
                case 2
                    res = min(min(F0));
                case 3
                    res = -max(max(abs(F0)));
                case 4
                    res = F0(1)-norm(F0(2:end));
                case 5
                    res = 2*F0(1)*F0(2)-norm(F0(3:end))^2
            end
        end
        data{j,1} = ['#' num2str(j)];
        data{j,2} = F.clauses{j}.symbolic;
        data{j,3} = lmiinfo{F.clauses{j}.type};
        data{j,4} = res;
        data{j,5} = F.clauses{j}.handle;
    end
end

res = [data{:,4}];