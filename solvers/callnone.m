function output = callnone(interfacedata)

% Author Johan Löfberg 
% $Id: callnone.m,v 1.2 2005-05-07 13:53:20 joloef Exp $

% Standard interface 
output.Primal      = zeros(length(interfacedata.c),1);
output.Dual        = [];
output.Slack       = [];
output.problem     = -2;
output.infostr     = yalmiperror(-2,'YALMIP');
output.solverinput = [];
output.solveroutput= [];
output.solvertime  = 0;