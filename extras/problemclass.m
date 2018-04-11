function p = problemclass(self)

% LP = 1;
% QP = 2;
% SOCP = 3;
% SDP = 4;
% MILP = 5;
% MIQP = 6;
% 
% possible = ones(1,6);
% 
% if nnz(model.variabletype)>0
%     possible([LP QP SOCP SDP MILP MIQP 
%  
model = self.problemclass;

combinatorial = (model.constraint.integer | model.constraint.binary | model.constraint.semicont | model.constraint.sos1 | model.constraint.sos2);

p = 'Unclassified';

if nnz(self.variabletype)==0
    % can be LP,QP,SDP,SOCP,MILP,MIQP
    if model.constraint.inequalities.semidefinite.linear
        p = 'SDP';
    elseif model.constraint.inequalities.secondordercone.linear | model.constraint.inequalities.rotatedsecondordercone.linear
        p = 'SOCP';
    else
        if model.objective.quadratic.convex == 1
            p = 'Convex QP';
        elseif model.objective.quadratic.nonconvex == 1
            p = 'Nonconvex QP';
        else
            p = 'LP';
        end
    end
    if combinatorial
        p = ['Mixed Integer ' p];
    end
end
