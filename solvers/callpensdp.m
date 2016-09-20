function output = callpensdpm(interfacedata);

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
ub      = interfacedata.ub;
lb      = interfacedata.lb;

% Bounded variables converted to constraints
if ~isempty(ub)
    [F_struc,K] = addStructureBounds(F_struc,K,ub,lb);
end

% *******************************
% CONVERT FROM INTERNAL FORMAT
% *******************************
penstruct = sedumi2pen(F_struc,K,c,x0);

% **************************
% COPY OPTIONS
% **************************
ops = struct2cell(options.pensdp);ops = [ops{1:end}];
penstruct.ioptions = ops(1:8);
penstruct.foptions = ops(9:end);
penstruct.ioptions(4) = max(0,min(3,options.verbose+1));
if penstruct.ioptions(4)==1
    penstruct.ioptions(4)=0;
end

% FIX
if penstruct.mconstr == 0
    penstruct.msizes = [];
end

if options.savedebug
    save pensdpdebug penstruct
end

%**************************
% CALL PENSDP
%**************************
showprogress('Calling PENSDP',options.showprogress);
solvertime = tic;
[x, fx, u, iresults, fresults, iflag] = pen(penstruct);
solvertime = toc(solvertime);

% Get dual variable (this must be possible to do easier...)
u = u(:);
D_struc = u(1:1:K.l);
if length(K.s)>0
    if K.s(1)>0
        pos = K.l+1;
        for i = 1:length(K.s)
            temp = zeros(K.s(i),K.s(i));
            vecZ = u(pos:pos+0.5*K.s(i)*(K.s(i)+1)-1);
            top = 1;
            for j = 1:K.s(i)
                len = K.s(i)-j+1;
                temp(j:end,j)=vecZ(top:top+len-1);top=top+len;
            end
            temp = (temp+temp');j = find(speye(K.s(i)));temp(j)=temp(j)/2;
            D_struc = [D_struc;temp(:)];
            pos = pos + (K.s(i)+1)*K.s(i)/2;
        end
    end
end

switch iflag 
case 0
    problem = 0;
case 1
    problem = 4;
case 2
    problem = 1;
case 3
    problem = 4;
case 4
    problem = 3;
case 5
    problem = 7;
case 6
    problem = 11;
case 7
    problem = 9;
otherwise
    problem = -1;
end    
infostr = yalmiperror(problem,'PENSDP/TOMLAB');	

if options.savesolveroutput
    solveroutput.f = f;
    solveroutput.x = x;
    solveroutput.u = u;
    solveroutput.iflag = iflag;
    solveroutput.niter = niter;
    solveroutput.feas = feas;
else
    solveroutput = [];
end
if options.savesolverinput
    solverinput.penstruct = penstruct;
else
    solverinput = [];
end

% Standard interface 
output = createOutputStructure(x(:),D_struc,[],problem,infostr,solverinput,solveroutput,solvertime);



