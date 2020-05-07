function output = callscs(model)
originalModel = model; 
% *********************************************
% Bounded variables converted to constraints
% N.B. Only happens when caller is BNB
% *********************************************
if ~isempty(model.ub)
    [model.F_struc,model.K] = addStructureBounds(model.F_struc,model.K,model.ub,model.lb);
end

% Write nonlinear functions using exponential cone operators, if possible
[model,output] = normalizeExponentialCone(model);
if output.problem
    return
end

% Map to SCS format
data.A = -model.F_struc(:,2:end);
data.b = full(model.F_struc(:,1));
data.c = full(model.c);
cones = [];
cones.f = model.K.f;
cones.l = model.K.l;
cones.q = model.K.q;
cones.s = model.K.s;
cones.ep = model.K.e;

params = model.options.scs;
params.verbose = model.options.verbose;

if params.eliminateequalities
    tempmodel.F_struc = [data.b -data.A];
    tempmodel.c = data.c;
    tempmodel.K = cones;
    [tempmodel,xsol,H] =  removeequalitiesinmodel(tempmodel);
    data.b = full(tempmodel.F_struc(:,1));
    data.A = -tempmodel.F_struc(:,2:end);
    data.c = full(tempmodel.c);
    cones = tempmodel.K;
else
    H = 1;
    xsol = 0;
end

% Extract lower diagonal form for new SCS format
if ~isempty(cones.s) && any(cones.s)
    sdpA = data.A(1+cones.l + cones.f+sum(cones.q):end,:);
    sdpb = data.b(1+cones.l + cones.f+sum(cones.q):end,:);
    expA = data.A(end-3*cones.ep+1:end,:);
    expb = data.b(end-3*cones.ep+1:end,:);
    data.A = data.A(1:cones.l + cones.f+sum(cones.q),:);    
    data.b = data.b(1:cones.l + cones.f+sum(cones.q),:);
    top = 1;
    for i = 1:length(cones.s)
        A = sdpA(top:top + cones.s(i)^2-1,:);
        b = sdpb(top:top + cones.s(i)^2-1,:);
        n = cones.s(i);
        ind = find(speye(n));
        b(ind) = b(ind)/sqrt(2);
        A(ind,:) = A(ind,:)/sqrt(2);
        ind = find(tril(ones(n)));
        A = A(ind,:);
        b = b(ind);
        data.A = [data.A;A];
        data.b = [data.b;b];
        top = top  + cones.s(i)^2;
    end
    data.A = [data.A;expA];
    data.b = [data.b;expb];
end

if model.options.savedebug
    save scsdebug data cones params
end

if model.options.showprogress;showprogress(['Calling ' model.solver.tag],model.options.showprogress);end
t = tic;
problem = 0;  
switch  model.solver.tag
    case 'scs-direct'
         [x_s,y_s,s,info] = scs_direct(data,cones,params);
    otherwise
        if params.gpu == true
            [x_s,y_s,s,info] = scs_gpu(data,cones,params);
        else
            [x_s,y_s,s,info] = scs_indirect(data,cones,params);
        end
end
solvertime = toc(t);

% Internal format. Only recover the original variables
Primal = H*x_s+xsol;
Primal = Primal(1:length(originalModel.c));
if ~isempty(model.evalMap)
    % No support for duals when exponential cones yet
    Dual = [];
else
    % Map to full format from tril
    Dual = y_s(1:cones.f+cones.l+sum(cones.q));
    if ~isempty(cones.s) && any(cones.s)        
        top = 1 + cones.f + cones.l + sum(cones.q);
        for i = 1:length(cones.s)
            n = cones.s(i);
            sdpdual = y_s(top:top + n*(n+1)/2-1,:);
            Z = zeros(n);
            Z(find(tril(ones(n)))) = sdpdual;
            Z = (Z + Z')/2;
            ind = find(speye(n));
            Z(ind) = Z(ind)/sqrt(2);
            Dual = [Dual;Z(:)];
            top = top  + n*(n+1)/2;
        end
    end
end

if nnz(data.c)==0 && isequal(info.status,'Unbounded/Inaccurate')
    info.status = 'Infeasible';
end

switch info.status
    case 'Solved'
        problem = 0;
    case 'Infeasible'
        problem = 1;
     case 'Unbounded'
        problem = 2;    
    otherwise
        status = 9;
end

infostr = yalmiperror(problem,model.solver.tag);

% Save ALL data sent to solver
if model.options.savesolverinput
    solverinput.data = data;
    solverinput.cones = cones;
    solverinput.params = params;  
else
    solverinput = [];
end

% Save ALL data from the solution?
if model.options.savesolveroutput
    solveroutput.x = x_s;
    solveroutput.y = y_s;
    solveroutput.s = s;
    solveroutput.info = info;
else
    solveroutput = [];
end

% Standard interface 
output = createOutputStructure(Primal,Dual,[],problem,infostr,solverinput,solveroutput,solvertime);

function [model,x0,H] = removeequalitiesinmodel(model)

if model.K.f > 0
    % Extract the inequalities
    A_equ = model.F_struc(1:model.K.f,2:end);
    b_equ = -model.F_struc(1:model.K.f,1);
    
    if 1
        % Find a basis for the column space of A_equ
        [L,U,P] = lu(A_equ');
        r = colspaces(L');
        AA = L';
        H1 = AA(:,r);
        H2 = AA(:,setdiff(1:size(AA,2),r));
        try
            x0 = P'*linsolve(full(L'),linsolve(full(U'),full(b_equ),struct('LT',1==1)),struct('UT',1==1));
        catch
            x0 = A_equ\b_equ;
        end
        % FIX : use L and U stupid!
        H = P'*[-H1\H2;eye(size(H2,2))];
        
    else
        H = null(full(A_equ));
        x0 = A_equ\b_equ;
    end
    
    model.c = H'*model.c;
    model.F_struc(:,1) = model.F_struc(:,1) + model.F_struc(:,2:end)*x0;
    model.F_struc = [model.F_struc(:,1) model.F_struc(:,2:end)*H];
    model.F_struc(1:model.K.f,:)=[];
    model.K.f = 0;
else
    H = speye(length(model.c));
    x0 = zeros(length(model.c),1);
end


function  [indx]=colspaces(A)
indx = [];
for i = 1:size(A,2)
    s = max(find(A(:,i)));
    indx = [indx s];
end
indx = unique(indx);