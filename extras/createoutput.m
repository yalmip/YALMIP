function output = createoutput(primal,dual,slack,problem,tag,input,solveroutput,time)

if nargin == 1
    problem = primal;
    primal = [];
    dual = [];
    slack = [];
    input = [];
    solveroutput = [];
    time = [];
end

output.Primal      = primal;
output.Dual        = dual;
output.Slack       = slack;
output.problem     = problem;
output.infostr     = yalmiperror(problem,tag);
output.solverinput = input;
output.solveroutput= solveroutput;
output.solvertime  = time;

