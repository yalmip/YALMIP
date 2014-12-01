function output = callscs(model)

% Retrieve needed data
options = model.options;
F_struc = model.F_struc;
c       = model.c;
K       = model.K;
ub      = model.ub;
lb      = model.lb;

% *********************************************
% Bounded variables converted to constraints
% N.B. Only happens when caller is BNB
% *********************************************
if ~isempty(ub)
    [F_struc,K] = addbounds(F_struc,K,ub,lb);
end

data.A = -F_struc(:,2:end);
data.b = full(F_struc(:,1));
data.c = full(c);
cones = K;

% Now add exponential cone information
if ~isempty(model.evalMap)
    % First check that we have only exponentials
    convexFunctions = [];
    concaveFunctions = [];
    vectorentropy = []; 
    vectorkullback = [];  
    vectorlogsumexp = [];
    for i = 1:length(model.evalMap)
        switch model.evalMap{i}.fcn
            case {'exp','pexp'}
                convexFunctions = [convexFunctions model.evalMap{i}.computes];               
            case {'log','plog','slog'}
                concaveFunctions = [concaveFunctions model.evalMap{i}.computes];
            case 'entropy'                              
                concaveFunctions = [concaveFunctions model.evalMap{i}.computes];
                if length(model.evalMap{i}.variableIndex) > 1
                    vectorentropy = [vectorentropy i];
                end
            case 'kullbackleibler'
                convexFunctions = [convexFunctions model.evalMap{i}.computes];
                if length(model.evalMap{i}.variableIndex) > 2
                    vectorkullback = [vectorkullback i];
                end
            case 'logsumexp'
                convexFunctions = [convexFunctions model.evalMap{i}.computes];                
                vectorlogsumexp = [vectorlogsumexp i];                   
            otherwise
                % Standard interface, return solver not applicable
                output = createoutput([],[],[],-4,model.solver.tag,[],[],0);
                return
        end
    end
    % Check that all exp/log enter in a convex fashion
    if model.K.f > 0
       if nnz(data.A(1:K.f,convexFunctions))>0 || nnz(data.A(1:K.f,concaveFunctions))>0
          output = createoutput([],[],[],-4,model.solver.tag,[],[],0);
             return
        end 
    end
    if any(data.c(convexFunctions) < 0) || any(data.c(concaveFunctions) > 0)
        output = createoutput([],[],[],-4,model.solver.tag,[],[],0);
        return
    end
    if model.K.l > 0
        if nnz(data.A(1+model.K.f:model.K.f+model.K.l,convexFunctions)<0) || nnz(data.A(1+model.K.f:model.K.f+model.K.l,concaveFunctions)>0)
             output = createoutput([],[],[],-4,model.solver.tag,[],[],0);
             return
        end
    end
    if sum(model.K.q) + sum(model.K.s) > 0
         if nnz(data.A(1+model.K.f+K.l:end,convexFunctions))>0 || nnz(data.A(1+model.K.f+K.l:end,concaveFunctions))>0
             output = createoutput([],[],[],-4,model.solver.tag,[],[],0);
             return
        end
    end
    % Scalarize entropy operator
    if ~isempty(vectorentropy)
        for i = vectorentropy(:)'
            involved = model.evalMap{i}.variableIndex;
            computes = model.evalMap{i}.computes;
            n = length(data.c);
            m = length(involved);
            data.A = [sparse(ones(1,m+1),[computes n+(1:m)],[-1 ones(1,m)],1,n+m);
                      data.A spalloc(size(data.A,1),m,0)];
            data.b = [0;data.b];      
            data.c = [data.c;zeros(m,1)];            
            cones.f = cones.f + 1;
            % The original strucuture now just relates to the first new term
            model.evalMap{i}.computes = n+1;
            model.evalMap{i}.variableIndex = involved(1);
            % and then we add some
            for j = 2:length(involved)
                model.evalMap{end+1} = model.evalMap{i};
                model.evalMap{end}.variableIndex = involved(j);
                model.evalMap{end}.computes = n+j;
            end
        end
    end

    % Scalarize kullbackleibler operator
    if ~isempty(vectorkullback)
        for i = vectorkullback(:)'
            involved = model.evalMap{i}.variableIndex;
            computes = model.evalMap{i}.computes;
            n = length(data.c);
            m = length(involved)/2;
            data.A = [sparse(ones(1,m+1),[computes n+(1:m)],[-1 ones(1,m)],1,n+m);
                      data.A spalloc(size(data.A,1),m,0)];
            data.b = [0;data.b];      
            data.c = [data.c;zeros(m,1)];            
            cones.f = cones.f + 1;
            % The original strucuture now just relates to the first new term
            model.evalMap{i}.computes = n+1;
            model.evalMap{i}.variableIndex = involved([1 m+1]);
            % and then we add some
            for j = 2:length(involved)/2
                model.evalMap{end+1} = model.evalMap{i};
                model.evalMap{end}.variableIndex = involved([j j+m]);
                model.evalMap{end}.computes = n+j;
            end
        end
    end
    % Scalarize logsumexp operator 
    % write log(expx1+expx2..)<= z as exp(x1-z)+exp(x2-z)+...<=1
    if ~isempty(vectorlogsumexp)
        for i = vectorlogsumexp(:)'
            involved = model.evalMap{i}.variableIndex;
            computes = model.evalMap{i}.computes;
            % Orignal #variables
            n = length(data.c);
            % Number of terms in logsumexp
            m = length(involved);           
            
            % exp(x1-z)+exp(x2-z)+...<=1
            data.A = [data.A(1:cones.f,:) spalloc(cones.f,m,0);
                      spalloc(1,n,0) ones(1,m);
                      data.A(cones.f+1:end,:) spalloc(cones.l+sum(cones.q)+sum(cones.s),m,0)];
            data.b = [data.b(1:cones.f);
                      1;
                      data.b(cones.f+1:end)];
            data.c = [data.c;zeros(m,1)];
            cones.l = cones.l + 1;   
            
            % The original strucuture now just relates to the first new term
            model.evalMap{i}.computes = n+1;
            model.evalMap{i}.variableIndex = [involved(1) computes];
            model.evalMap{i}.fcn = 'expdiff';
            % and then we add some new
            for j = 2:length(involved)
                model.evalMap{end+1} = model.evalMap{i};
                model.evalMap{end}.variableIndex = [involved(j) computes];
                model.evalMap{end}.computes = n+j;
            end
        end
    end
    
    % Describe all exponential cones
    m = length(model.evalMap);        
    cones.ep = m;
    dataAi = [];
    dataAj = [];
    dataAk = [];
    datab = [];
    for i = 1:m
        switch model.evalMap{i}.fcn
            case 'exp'
                % 1*exp(xv/1) <= xc
                % y*exp(x/y)  <= z.  y new variable, xv the original variable, and xc the "computed"
                x = model.evalMap{i}.variableIndex;
                z = model.evalMap{i}.computes;
                data.A = [data.A;sparse([1;3],[x z],-1,3,size(data.A,2))];
                data.b = [data.b;[0;1;0]];
            case 'pexp'
                % xv(1)*exp(xv(2)/xv(1)) <= xc    
                x = model.evalMap{i}.variableIndex(2);
                y = model.evalMap{i}.variableIndex(1);
                z = model.evalMap{i}.computes;
                data.A = [data.A;sparse([1;2;3],[x y z],[-1 -1 -1],3,size(data.A,2))];
                data.b = [data.b;[0;0;0]]; 
            case 'log'
                % log(xv) >= xc i.e. xv >= exp(xc/1)*1
                z = model.evalMap{i}.variableIndex;
                x = model.evalMap{i}.computes;
                data.A = [data.A;sparse([1;3],[x z],-1,3,size(data.A,2))];
                data.b = [data.b;[0;1;0]];
            case 'plog'
                % xv(1)log(xv(2)/xv(1))>=xc i.e.
                % -xc >= -xv(1)log(xv(2)/xv(1))
                z = model.evalMap{i}.computes;
                y = model.evalMap{i}.variableIndex(1);
                x = model.evalMap{i}.variableIndex(2);
                data.A = [data.A;sparse([1;2;3],[x y z],[1 -1 1],3,size(data.A,2))];
                data.b = [data.b;[0;0;0]];
            case 'slog'
                % log(1+xv) >= xc i.e. (1+xv) >= exp(xc/1)*1
                z = model.evalMap{i}.variableIndex;
                x = model.evalMap{i}.computes;
                data.A = [data.A;sparse([1;3],[x z],-1,3,size(data.A,2))];
                data.b = [data.b;0;1;1];
            case 'entropy'
                % -xv*log(xv)>=xc i.e. 1 >= exp(xc/xv)*xv
                x = model.evalMap{i}.computes;
                y = model.evalMap{i}.variableIndex;
                data.A = [data.A;sparse([1;2],[x y],[-1 -1],3,size(data.A,2))];
                data.b = [data.b;[0;0;1]];          
            case 'kullbackleibler'
                % -xv1*log(xv1/xv2)>=-xc i.e. xv2 >= exp(-xc/xv1)*xv1
                x = model.evalMap{i}.computes;
                y = model.evalMap{i}.variableIndex(1);
                z = model.evalMap{i}.variableIndex(2);
                data.A = [data.A;sparse([1;2;3],[x y z],[1 -1 -1],3,size(data.A,2))];
                data.b = [data.b;[0;0;0]];
            case 'expdiff'
                % exp(xv(1)-xv(2)) <= xc
                x = model.evalMap{i}.variableIndex;
                z = model.evalMap{i}.computes;
                data.A = [data.A;sparse([1;1;3],[x(1) x(2) z],[-1 1 -1],3,size(data.A,2))];
                data.b = [data.b;[0;1;0]];                                       
            otherwise
        end
    end
    if ~isempty(dataAi)
        data.A = [data.A;sparse(dataAi,dataAj,dataAk,3*length(dataAi)/3,size(data.A,2))];
        data.b = [data.b;datab];
    end            
end

if options.showprogress;showprogress(['Calling ' model.solver.tag],options.showprogress);end

params = options.scs;
params.verbose = options.verbose;

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

if options.savedebug
    save scsdebug data cones
end

tic
problem = 0;  
switch  model.solver.tag
    case 'scs-direct'
         [x_s,y_s,s,info] = scs_direct(data,cones,params);
    otherwise
        [x_s,y_s,s,info] = scs_indirect(data,cones,params);
end
solvertime = toc;

% Internal format. Only recover the original variables
Primal = H*x_s+xsol; 
Primal = Primal(1:length(model.c));
if isempty(model.evalMap)
    Dual = [];
else
    Dual = y_s;
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
if options.savesolverinput
    solverinput.data = data;
    solverinput.cones = cones;
    solverinput.params = params;  
else
    solverinput = [];
end

% Save ALL data from the solution?
if options.savesolveroutput
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