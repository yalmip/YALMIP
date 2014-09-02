function solver = findapplicablesolvers(options,ProblemClass,solvers,socp_are_really_qc);
%findapplicablesolvers Internal function to select solver based on problem category

% ************************************************
% Prune based on objective
% ************************************************
if ProblemClass.objective.sigmonial
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)
        keep(i) = solvers(i).objective.sigmonial;                         
    end
    solvers = solvers(find(keep));
end    
if ProblemClass.objective.polynomial
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)
        keep(i) =  solvers(i).constraint.equalities.quadratic | solvers(i).constraint.inequalities.elementwise.quadratic.nonconvex | solvers(i).objective.polynomial | solvers(i).objective.sigmonial;            
    end
    solvers = solvers(find(keep));
end  
if ProblemClass.objective.quadratic.nonconvex
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)
        keep(i) = solvers(i).objective.polynomial | solvers(i).objective.sigmonial | solvers(i).objective.quadratic.nonconvex;        
    end
    solvers = solvers(find(keep));
end  
if ProblemClass.objective.quadratic.convex
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)
        direct = solvers(i).objective.polynomial | solvers(i).objective.sigmonial | solvers(i).objective.quadratic.nonconvex | solvers(i).objective.quadratic.convex;
        indirect = solvers(i).constraint.inequalities.semidefinite.linear | solvers(i).constraint.inequalities.secondordercone;
        if direct | indirect
            keep(i)=1;
        else
            keep(i)=0;
        end
    end
    solvers = solvers(find(keep));
end  
if ProblemClass.objective.linear
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)
        keep(i) = solvers(i).objective.polynomial | solvers(i).objective.sigmonial | solvers(i).objective.quadratic.nonconvex | solvers(i).objective.quadratic.convex | solvers(i).objective.linear;
    end
    solvers = solvers(find(keep));
end  


% ******************************************************
% Prune based on rank constraints
% ******************************************************
if ProblemClass.constraint.inequalities.rank
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.rank;
    end
    solvers = solvers(find(keep));
end  

% ******************************************************
% Prune based on semidefinite constraints
% ******************************************************
if ProblemClass.constraint.inequalities.semidefinite.sigmonial
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.semidefinite.sigmonial;
    end
    solvers = solvers(find(keep));
end  
if ProblemClass.constraint.inequalities.semidefinite.polynomial
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.semidefinite.sigmonial |  solvers(i).constraint.inequalities.semidefinite.polynomial;
    end
    solvers = solvers(find(keep));
end  
if ProblemClass.constraint.inequalities.semidefinite.quadratic
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.semidefinite.sigmonial |  solvers(i).constraint.inequalities.semidefinite.polynomial | solvers(i).constraint.inequalities.semidefinite.quadratic;
    end
    solvers = solvers(find(keep));
end  
if ProblemClass.constraint.inequalities.semidefinite.linear
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.semidefinite.sigmonial |  solvers(i).constraint.inequalities.semidefinite.polynomial | solvers(i).constraint.inequalities.semidefinite.quadratic | solvers(i).constraint.inequalities.semidefinite.linear;
    end        
    solvers = solvers(find(keep));
end  

% ******************************************************
% Prune based on second order cone constraints
% ******************************************************
if ProblemClass.constraint.inequalities.secondordercone & ~socp_are_really_qc
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
         keep(i) = solvers(i).constraint.inequalities.secondordercone | solvers(i).constraint.inequalities.semidefinite.linear;
    end
    solvers = solvers(find(keep));
end  
if ProblemClass.constraint.inequalities.rotatedsecondordercone
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.rotatedsecondordercone | solvers(i).constraint.inequalities.secondordercone | solvers(i).constraint.inequalities.semidefinite.linear;
    end
    solvers = solvers(find(keep));
end  

% ******************************************************
% Prune based on element-wise inequality constraints
% ******************************************************
if ProblemClass.constraint.inequalities.elementwise.sigmonial
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.elementwise.sigmonial;            
    end
    solvers = solvers(find(keep));
end  
if ProblemClass.constraint.inequalities.elementwise.polynomial
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.elementwise.quadratic.nonconvex | solvers(i).constraint.inequalities.elementwise.sigmonial |  solvers(i).constraint.inequalities.elementwise.polynomial;
    end
    solvers = solvers(find(keep));
end  
if ProblemClass.constraint.inequalities.elementwise.quadratic.nonconvex
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.elementwise.sigmonial |  solvers(i).constraint.inequalities.elementwise.polynomial | solvers(i).constraint.inequalities.elementwise.quadratic.nonconvex;
    end
    solvers = solvers(find(keep));
end 
if ProblemClass.constraint.inequalities.elementwise.quadratic.convex | (ProblemClass.constraint.inequalities.secondordercone & socp_are_really_qc)
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.elementwise.sigmonial |  solvers(i).constraint.inequalities.elementwise.polynomial | solvers(i).constraint.inequalities.elementwise.quadratic.nonconvex | solvers(i).constraint.inequalities.elementwise.quadratic.convex | solvers(i).constraint.inequalities.secondordercone | solvers(i).constraint.inequalities.semidefinite.linear;
    end
    solvers = solvers(find(keep));
end 
if ProblemClass.constraint.inequalities.elementwise.linear
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.elementwise.sigmonial |  solvers(i).constraint.inequalities.elementwise.polynomial | solvers(i).constraint.inequalities.semidefinite.quadratic | solvers(i).constraint.inequalities.elementwise.linear;
    end
    solvers = solvers(find(keep));
end  

% ******************************************************
% Prune based on element-wise constraints
% ******************************************************
if ProblemClass.constraint.equalities.sigmonial
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.inequalities.elementwise.sigmonial | solvers(i).constraint.equalities.sigmonial;            
    end
    solvers = solvers(find(keep));
end  
if ProblemClass.constraint.equalities.polynomial
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)    
        indirect = solvers(i).constraint.inequalities.elementwise.quadratic.nonconvex | solvers(i).constraint.inequalities.elementwise.sigmonial |  solvers(i).constraint.inequalities.elementwise.polynomial;
        indirect = indirect | solvers(i).constraint.inequalities.elementwise.sigmonial | solvers(i).constraint.inequalities.elementwise.polynomial;
        direct = solvers(i).constraint.equalities.sigmonial |  solvers(i).constraint.equalities.polynomial;
        keep(i) = direct | indirect;
    end
    solvers = solvers(find(keep));
end  
if ProblemClass.constraint.equalities.quadratic
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        indirect = solvers(i).constraint.inequalities.elementwise.sigmonial | solvers(i).constraint.inequalities.elementwise.polynomial | solvers(i).constraint.inequalities.elementwise.quadratic.nonconvex;
        direct = solvers(i).constraint.equalities.sigmonial |  solvers(i).constraint.equalities.polynomial |  solvers(i).constraint.equalities.quadratic;
        keep(i) = direct | indirect;
    end
    solvers = solvers(find(keep));
end 
if ProblemClass.constraint.equalities.linear
    keep = ones(length(solvers),1);
    for i = 1:length(solvers) 
        indirect = solvers(i).constraint.inequalities.elementwise.linear | solvers(i).constraint.inequalities.elementwise.sigmonial |  solvers(i).constraint.inequalities.elementwise.polynomial;
        direct = solvers(i).constraint.equalities.linear | solvers(i).constraint.equalities.sigmonial |  solvers(i).constraint.equalities.polynomial;
        keep(i) = direct | indirect;
    end
    solvers = solvers(find(keep));
end  


% ******************************************************
% Discrete data
% ******************************************************
if ProblemClass.constraint.integer
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
        keep(i) = solvers(i).constraint.integer;
    end
    solvers = solvers(find(keep));
end  
if ProblemClass.constraint.binary
    keep = ones(length(solvers),1);
    for i = 1:length(solvers)                      
         keep(i) = solvers(i).constraint.integer | solvers(i).constraint.binary;            
    end
    solvers = solvers(find(keep));
end  

if isempty(solvers)
    solver = [];
else
    solver=solvers;
 %   solver = solvers(1);
end

% if isempty(solver)
%     if length(options.solver)>0 % User selected available solver, but it is not applicable
%         problem = -4;
%     else
%         problem = -2;
%     end
% end


% FIX : Hack when chosing the wrong fmincon thingy
if ~isempty(solver)
    if (length(options.solver)==0 | isequal(options.solver,'fmincon')) & isequal(solver.tag,'fmincon') & isequal(solver.version,'geometric')
        if ~(ProblemClass.objective.sigmonial | ProblemClass.constraint.inequalities.elementwise.sigmonial)
            solver.version = 'standard';
            solver.call    = 'callfmincon';
            solver.objective.linear = 1;
            solver.objective.quadratic.convex = 1;
            solver.objective.quadratic.nonconvex = 1;
            solver.objective.polynomial = 1;
            solver.objective.sigmonial = 1;
            solver.constraint.equalities.elementwise.linear = 1;
            solver.constraint.equalities.elementwise.quadratic.convex = 1;
            solver.constraint.equalities.elementwise.quadratic.nonconvex = 1;
            solver.constraint.equalities.elementwise.polynomial = 1;
            solver.constraint.equalities.elementwise.sigmonial = 1;
            solver.constraint.inequalities.elementwise.linear = 1;
            solver.constraint.inequalities.elementwise.quadratic.convex = 1;
            solver.constraint.inequalities.elementwise.quadratic.nonconvex = 1;
            solver.constraint.inequalities.elementwise.polynomial = 1;
            solver.constraint.inequalities.elementwise.sigmonial = 1;
            solver.constraint.inequalities.semidefinite.linear = 1;
            solver.constraint.inequalities.semidefinite.quadratic = 1;
            solver.constraint.inequalities.semidefinite.polynomial = 1;
            solver.dual = 1;
        end
    end
end


