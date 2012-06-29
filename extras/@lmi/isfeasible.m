function [feasible,feaslistLMI] = isfeasible(F,tol)

% Author Johan Löfberg 
% $Id: isfeasible.m,v 1.4 2005-12-17 12:55:31 joloef Exp $   

if nargin == 1
    tol = 0;
end
feaslistLMI = checkset(F);
feasible = all(feaslistLMI >= -tol);

return

% Check if solution avaliable
currsol = evalin('caller','yalmip(''getsolution'')');
if isempty(currsol)
    feasible = 0;
end

nlmi = size(F.clauses,2);
if (nlmi == 0)
    feasible = 1;
    feaslistLMI = [];
    return
end

feaslistLMI = zeros(nlmi,1);
if nlmi>0
    for j = 1:nlmi
        F0 = double(F.clauses{j}.data);
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
        feaslistLMI(j) = res;
    end
end

if nargin==1
    if (isempty(feaslistLMI) | min(feaslistLMI)>=0)
        feasible = 1;
    else
        feasible = 0;
    end
end
