function val = value(P)

val = 0;
for i = 1:length(P.Constraint)
    risk_i = 1 - value(P.Risk{i});
    val = val + value(P.Weight{i})*(P.Offset{i}+risk_i);
end
