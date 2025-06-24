function sys = display(P)
if P.joint
    display(['Joint probability expression (' num2str(length(P.Constraint)) ' term)'])
else
    display(['Probability expression (' num2str(length(P.Constraint)) ' term)'])
end