function output = callsdpagmp(interfacedata)

% Path to executable
path2sdpagmp = '/usr/local/bin/';

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
ub      = interfacedata.ub;
lb      = interfacedata.lb;

% Check if param.sdpa file exists in pwd
cleanup = 0;
if ~exist([pwd,filesep,'param.sdpa'],'file')
    cleanup = 1;
    % write it with default parameters, otherwise failure!
    fID = fopen([pwd,filesep,'param.sdpa'],'w');
    fprintf(fID,'200         unsigned int    maxIteration;                 \n');
    fprintf(fID,'1.0E-25     double          0.0 < epsilonStar;            \n');
    fprintf(fID,'1.0E6       double          0.0 < lambdaStar;             \n');
    fprintf(fID,'2.0         double          1.0 < omegaStar;              \n');
    fprintf(fID,'-1.0E25     double          lowerBound;                   \n');
    fprintf(fID,'1.0E25      double          upperBound;                   \n');
    fprintf(fID,'0.1         double          0.0 <= betaStar <  1.0;       \n');
    fprintf(fID,'0.2         double          0.0 <= betaBar  <  1.0, betaStar <= betaBar;\n');
    fprintf(fID,'0.7         double          0.0 < gammaStar <  1.0;       \n');
    fprintf(fID,'1.0E-25     double          0.0 < epsilonDash;            \n');
    fprintf(fID,'200         precision;                                    \n');
    fclose(fID);
end

% Bounded variables converted to constraints
if ~isempty(ub)
    % addbounds was renamed somewhere along the way in YALMIP?
    try
        [F_struc,K] = addbounds(F_struc,K,ub,lb);
    catch
        [F_struc,K] = addStructureBounds(F_struc,K,ub,lb);
    end
end

% Convert from internal (sedumi) format
[mDIM,nBLOCK,bLOCKsTRUCT,c,F] = sedumi2sdpa(F_struc,c,K);

if options.verbose==0
    options.sdpa.print = 'no';
else
    options.sdpa.print = 'display';
end

if options.savedebug
    ops = options.sdpa;
    save sdpadebug mDIM nBLOCK bLOCKsTRUCT c F ops
end

if options.showprogress;
    showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALL SDPA-GMP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export to SDPA-GMP, solve and inport results
FF = cellfun(@full,F,'UniformOutput',false);
[nmax,mmax] = size(FF);
BS = abs(bLOCKsTRUCT);
for nn=1:nmax               % loop to construc SDPA problem
    for mm=1:mmax
        if isempty(FF{nn,mm})
            FF{nn,mm} = zeros(BS(nn));
        end
    end
end

tic
inputSDPA  = ['sdpagmp_in.dat-s'];
outputSDPA = ['sdpagmp_out.out'];
header     = 'Input from YALMIP';

gensdpagmpfile(inputSDPA,mDIM,nBLOCK,bLOCKsTRUCT,c,FF,header);          % write SDPA-GMP input file
if options.verbose==0
    system(['echo ',repmat('+',1,100),' >> sdpagmp.log']);
    system([path2sdpagmp,'sdpa_gmp ',inputSDPA,' ',outputSDPA,' >> sdpagmp.log']);          % solve SDP
else
    system(['echo ',repmat('+',1,100)]);
    system([path2sdpagmp,'sdpa_gmp ',inputSDPA,' ',outputSDPA]);          % solve SDP
end
% import result
[objVal,x,X,Y,INFO] = sdpagmp_read_output(outputSDPA,full(mDIM),full(nBLOCK),full(bLOCKsTRUCT));
solvertime = toc;

% Clean up tmp files created in this directory
delete(inputSDPA);
delete(outputSDPA);
if cleanup
    delete('param.sdpa');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From here onwards, like in YALMIP native callsdpa
% Create variables in YALMIP internal format
Primal = x;

Dual = [];
for i = 1:length(Y)
    Dual = [Dual;Y{i}(:)];
end

Slack = [];
if options.saveduals
    for i = 1:length(X)
        Slack = [Slack;X{i}(:)];
    end
end

switch (INFO.phasevalue)
    case 'pdOPT'
        problem = 0;
    case {'noINFO','pFEAS','dFEAS'}
        problem = 3;
    case {'pdFEAS'}
        problem = 4;
    case 'pFEAS_dINF'
        problem = 2;
    case 'pINF_dFEAS'
        problem = 1;
    case 'pUNBD'
        problem = 2;
    case 'dUNBD'
        problem = 1;
    case 'pdINF'
        problem = 12;
    otherwise
        problem = -1;
end
infostr = yalmiperror(problem,interfacedata.solver.tag);

if options.savesolveroutput
    solveroutput.objVal = objVal;
    solveroutput.x = x;
    solveroutput.X = X;
    solveroutput.Y = Y;
    solveroutput.INFO = INFO;
else
    solveroutput = [];
end

if options.savesolverinput
    solverinput.mDIM = mDIM;
    solverinput.nBLOCK=nBLOCK;
    solverinput.bLOCKsTRUCT=bLOCKsTRUCT;
    solverinput.c=c;
    solverinput.F=F;
else
    solverinput = [];
end

% Standard interface
output = createOutputStructure(Primal,Dual,[],problem,infostr,solverinput,solveroutput,solvertime);

























