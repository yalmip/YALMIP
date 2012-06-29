function Fcont = logic2cont(Flogic)

% Author Johan Löfberg
% $Id: logic2cont.m,v 1.1 2004-08-06 14:23:37 johanl Exp $

Fcont = set([]);
M = 1e3;
% How many elements in this clause
n = size( Flogic.clauses{1},2);
% Are all simple letters?
all_letters = 1;
j = 1;
while all_letters & j<=n    
    if isa(Flogic.clauses{1}{j}.data,'sdpvar')
        all_letters = all_letters &  is(Flogic.clauses{1}{j}.data,'logic');
    else
        all_letters = 0;
    end
    j = j + 1;
end

if all_letters
    sum_or = 0;
    % Ah, simple constraints, just add them
    for j = 1:n
        sum_or = sum_or + Flogic.clauses{i}{j};
    end
else
    sum_or = 0;
    delta = binvar(size( Flogic.clauses{i},2),1);
    for j = 1:size( Flogic.clauses{i},2)
        clause_i_j = Flogic.clauses{i}{j};
        if isa(clause_i_j,'constraint')
            delta = binvar(1,1);
            F = F + set(sdpvar(clause_i_j)+M*(1-delta)>0);
            sum_or = sum_or + delta;
        elseif isa(clause_i_j,'sdpvar')
            F = F + set(clause_i_j);
            sum_or = sum_or + clause_i_j;
        end
    end

end
Fcont = set(sum_or>1);