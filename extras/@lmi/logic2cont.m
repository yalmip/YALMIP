function Fcont = logic2cont(Flogic)

Fcont = ([]);
M = 1e3;
% How many elements in this clause
Flogic = flatten(Flogic);
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
            F = F + (sdpvar(clause_i_j)+M*(1-delta)>=0);
            sum_or = sum_or + delta;
        elseif isa(clause_i_j,'sdpvar')
            F = F + (clause_i_j>=0);
            sum_or = sum_or + clause_i_j;
        end
    end

end
Fcont = (sum_or>=1);