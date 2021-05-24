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

scsmodel = yalmip2scs(model);

if model.options.scs.eliminateequalities
    tempmodel.F_struc = [scsmodel.data.b -scsmodel.data.A];
    tempmodel.c = scsmodel.data.c;
    tempmodel.K = scsmodel.cones;
    [tempmodel,xsol,H] =  removeequalitiesinmodel(tempmodel);
    scsmodel.data.b = full(tempmodel.F_struc(:,1));
    scsmodel.data.A = -tempmodel.F_struc(:,2:end);
    scsmodel.data.c = full(tempmodel.c);
    scsmodel.cones = tempmodel.K;
else
    H = 1;
    xsol = 0;
end


if model.options.savedebug
    save scsdebug scsmodel
end

if model.options.showprogress;showprogress(['Calling ' model.solver.tag],model.options.showprogress);end
t = tic;
problem = 0;  
switch  model.solver.tag
    case 'scs-direct'
         [x_s,y_s,s,info] = scs_direct(scsmodel.data,scsmodel.cones,scsmodel.param);
    otherwise
        if scsmodel.param.gpu == true
            [x_s,y_s,s,info] = scs_gpu(scsmodel.data,scsmodel.cones,scsmodel.param);
        else
            [x_s,y_s,s,info] = scs_indirect(scsmodel.data,scsmodel.cones,scsmodel.param);
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
    Dual = y_s(1:scsmodel.cones.f+scsmodel.cones.l+sum(scsmodel.cones.q));
    if ~isempty(scsmodel.cones.s) && any(scsmodel.cones.s)        
        top = 1 + scsmodel.cones.f + scsmodel.cones.l + sum(scsmodel.cones.q);
        for i = 1:length(scsmodel.cones.s)
            n = scsmodel.cones.s(i);
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

if nnz(scsmodel.data.c)==0 && isequal(info.status,'Unbounded/Inaccurate')
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

% Save ALL data sent to solver
if model.options.savesolverinput
    solverinput = scsmodel;
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
output = createOutputStructure(Primal,Dual,[],problem,model.solver.tag,solverinput,solveroutput,solvertime);

function [model,x0,H] = removeequalitiesinmodel(model)

if any(model.K.f)
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