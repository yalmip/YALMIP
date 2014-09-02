function sol = solveranksos(F,obj,options,ranks,BlockedQ)

Frank = ([]);
for i = 1:length(ranks)
    if ~isinf(ranks(i))
        Frank = Frank + (rank(BlockedQ{i}{1}) <= ranks(i));
    end
end
% rank adds the pos.def constraints again!!, so we remove them
check = ones(length(F),1);
keep  = ones(length(F),1);
for i = 1:length(BlockedQ)
    for j = 1:length(BlockedQ{i})
        Qij = BlockedQ{i}{j};
        for k = find(check)'
            if isequal(Qij,sdpvar(F(k)))
                keep(k)  = 0;
                check(k) = 0;
            end
        end
    end
end
% Let's hope LMIRANK is there
sol =  solvesdp(F(find(keep)) + Frank,[],sdpsettings(options,'solver',''));