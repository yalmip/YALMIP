function output = callspecsdp(interfacedata)

% Author Johan Löfberg 
% $Id: callspecsdp.m,v 1.3 2005-05-07 13:53:20 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
ub      = interfacedata.ub;
lb      = interfacedata.lb;
x0      = interfacedata.x0;

% *********************************************
% Bounded variables converted to constraints
% N.B. Only happens when caller is BNB
% *********************************************
if ~isempty(ub)
    [F_struc,K] = addbounds(F_struc,K,ub,lb);
end

% *********************************************
% Lower triangular parts on semi-definite parts
% *********************************************
trilindicies = [];
top = 0;
Dims = [];
Nlmi = 0;
nbvar = length(c);
if K.l>0
    trilindicies = (1:K.l)';
    top = top + K.l;
    Dims = [Dims ones(1,K.l)];
    Nlmi = Nlmi + K.l;
end
if K.s(1)>0
    for i = 1:length(K.s)
        trilindicies = [trilindicies;top + find(tril(ones(K.s(i))))];
        top = top + K.s(i)^2;
    end
    Nlmi = Nlmi + length(K.s);
    Dims = [Dims K.s];
end
A0 = full(-F_struc(trilindicies,1));
AA = -F_struc(trilindicies,2:end);

% YALMIP sends empty initials
if isempty(x0)
    x0=zeros(nbvar,1);
end

% Apkarians code to get indicies in sparse format
% Note, code to remove non-used variables is 
% skipped, YALMIP takes care of this (I think)
nrowA=size(AA,1);
[irowA,icolA,AAval]=find(AA);
nonzerA=length(icolA);
nbvar=max(icolA);
nnnn=[nrowA nonzerA nbvar Nlmi];

% *********************************************
% Get options
% *********************************************
ops     = struct2cell(options.specsdp);ops = [ops{1:end}];
opts    = ops(1:6);
penopts = ops(7:end);

% *********************************************
% Call Apkarians solver
% *********************************************
solvertime = clock; 
showprogress('Calling Apkarian',options.showprogress);
problem = 0;
[x,laug,lmax]=SPSDPLSCX(nnnn,A0,AAval,irowA,icolA,Dims,c,opts,penopts,x0);
solvertime = etime(clock,solvertime);

% *********************************************
% Dual variables not available
% *********************************************
D_struc = [];

% *********************************************
% Error codes not available
% *********************************************
problem = 13;
infostr = yalmiperror(problem,'SPECSDP');

% Save ALL data sent to solver
if options.savesolverinput
    solverinput = [];
else
    solverinput = [];
end

% Save ALL data from the solution?
if options.savesolveroutput
	solveroutput = [];
else
	solveroutput = [];
end

% Standard interface 
output.Primal      = x(:);
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;