function output = callnone(interfacedata)

% Standard interface 
output.Primal      = zeros(length(interfacedata.c),1);
output.Dual        = [];
output.Slack       = [];
output.problem     = -2;
output.infostr     = yalmiperror(-2,'YALMIP');
output.solverinput = [];
output.solveroutput= [];
output.solvertime  = 0;