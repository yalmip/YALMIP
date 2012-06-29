function [F,ProblemClass] = modelComplementarityConstraints(F,solver,ProblemClass)

i = find(is(F,'complementarity'));
Fc = F(i);
F(i)=[];
for i = 1:length(Fc)
    X = sdpvar(Fc(i));
    C1 = X(:,1);
    C2 = X(:,2);
    if solver.constraint.equalities.polynomial == 1
        [M,m]=derivebounds(C1);
        F = [F,C1>=0, C2>=0, C1'*C2 == 0];
        ProblemClass.constraint.equalities.polynomial = 1;
    elseif (solver.constraint.binary == 1) | (solver.constraint.integer == 1)
        [M1,m1]=derivebounds(C1);
        [M2,m2]=derivebounds(C2);
        delta = struct(Fc(i)).clauses{1}.extra.indicators;
        % delta = binvar(length(C1),1);        
        F = [F, C1>=0, C2>=0, C1<= M1.*delta, C2 <= M2.*(1-delta)];
        ProblemClass.constraint.binary = 1;
        ProblemClass.constraint.inequalities.elementwise.linear = 1;        
    end
end
