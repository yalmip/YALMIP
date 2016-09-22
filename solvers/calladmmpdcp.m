function output = callsedumi(model)

% Retrieve needed data
options = model.options;
F_struc = model.F_struc;
c       = model.c;
K       = model.K;
ub      = model.ub;
lb      = model.lb;

% Create the parameter structure
pars = options.sedumi;
pars.fid = double(options.verbose);

% *********************************************
% Bounded variables converted to constraints
% N.B. Only happens when caller is BNB
% *********************************************
if ~isempty(ub)
    [F_struc,K] = addStructureBounds(F_struc,K,ub,lb);
end

% Avoid bug (by design?) in SeDuMi on 1D second order cones
if any(K.q == 2)
    [F_struc,K] = fix1Dqcone(F_struc,K);
end

options.admmpdcp.verbose = options.verbose;
aK.f = K.f;
aK.l = K.l;
aK.q = K.q;
aK.s = K.s;
K = aK;

if options.savedebug
    save admmpdcp F_struc c K pars
end

% *********************************************
% Call admmPDCP
% *********************************************
if options.showprogress;showprogress(['Calling ' model.solver.tag],options.showprogress);end

problem = 0;  
solvertime = tic;
try       
    [x_s,y_s,info] = admmPDCP(-F_struc(:,2:end),-c,F_struc(:,1),K,options.admmpdcp);
    Primal = y_s;
    Dual   = x_s;
catch    
    [sol,info] = admmPDCP(-F_struc(:,2:end),-c,F_struc(:,1),K,options.admmpdcp);    
    Primal = [];  
    Dual = [];      
end
solvertime = toc(solvertime);

problem = 0;

infostr = yalmiperror(problem,model.solver.tag);

% Save ALL data sent to solver
if options.savesolverinput
    solverinput.A = -F_struc(:,2:end);
    solverinput.c = F_struc(:,1);
    solverinput.b = -c;
    solverinput.K = K;
    solverinput.ops = options.admmpdcp;
else
    solverinput = [];
end

% Save ALL data from the solution?
if options.savesolveroutput
    solveroutput.x = x_s;
    solveroutput.y = y_s;
    solveroutput.info = info;
else
    solveroutput = [];
end

% Standard interface 
output = createOutputStructure(Primal,Dual,[],problem,infostr,solverinput,solveroutput,solvertime);

 
function [new_F_struc,K] = fix1Dqcone(F_struc,K); 

new_F_struc = F_struc(1:(K.l+K.f),:);
top = K.f + K.l + 1;
for i = 1:length(K.q)
    if K.q(i) == 2
      new_F_struc = [new_F_struc;F_struc(top,:);F_struc(top+1,:);F_struc(top,:)*0];
      K.q(i) = 3;
      top = top + 2;
    else
      new_F_struc = [new_F_struc;F_struc(top:top+K.q(i)-1,:)];
      top = top + K.q(i);
    end
end
new_F_struc = [new_F_struc;F_struc(top:end,:)];