function output = createOutputStructure(Primal,Dual,Slack,problem,infostr,solverinput,solveroutput,solvertime)
if nargin == 1
    % Quick call for empty with just an error code
    output.problem = Primal;
    output.Primal      = [];
    output.Dual        = [];
    output.Slack       = [];
    output.infostr     = '';
    output.solverinput = [];
    output.solveroutput= [];
    output.solvertime  = 0;
    return
end
% Standard interface 
output.Primal      = Primal;
output.Dual        = Dual;
output.Slack       = [];
output.problem     = problem;
output.infostr     = yalmiperror(problem,infostr);
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;