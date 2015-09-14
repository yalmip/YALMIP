function output = callscipnl(model)

% This sets up everything and more. Can be simplified significantly since
% baron handles its own computational tree etc
model = yalmip2nonlinearsolver(model);

% [Anonlinear*f(x) <= b;Anonlinear*f(x) == b]
cu = full([model.bnonlinineq;model.bnonlineq]);
cl = full([repmat(-inf,length(model.bnonlinineq),1);model.bnonlineq]);
Anonlinear = [ model.Anonlinineq; model.Anonlineq];

% Create string representing objective
obj = createmodelstring(model.c,model);
obj = strrep(obj,'sqrtm_internal','sqrt');
if nnz(model.Q)>0
    obj = [obj '+' createQstring(model.Q,model)];
end
if model.f > 0
    obj = [obj '+' num2str(model.f)];
end
if isempty(obj)
    obj = [];
else
    if obj(1)=='+'
        obj = obj(2:end);
    end
% Append quadratic term
    obj = ['@(x) ' obj];
    obj = eval(obj);
end

% Create string representing nonlinear constraints
if length(cu)>0
    con = '[';
    remove = [];
    for i = 1:length(cu)
        if isinf(cl(i)) & isinf(cu(i))
            remove = [remove;i];
        else
            con = [con createmodelstring(Anonlinear(i,:),model) ';'];
        end
    end
    cl(remove) = [];
    cu(remove) = [];
    con = [con ']'];
    con = strrep(con,'sqrtm_internal','sqrt');
    con = strrep(con,'slog(x','log(1+x');
    con = ['@(x) ' con];
    con = eval(con);
else
    cu = [];
    cl = [];
    con = [];
end

% Linear constraints
ru = full([model.b;model.beq]);
rl = full([repmat(-inf,length(model.b),1);model.beq]);
A =  [model.A; model.Aeq];
if length(ru)>0
    remove = find(isinf(rl) & isinf(ru));
    A(remove,:)=[];
    rl(remove) = [];
    ru(remove) = [];
end

lb = model.lb;
ub = model.ub;

xtype = [];
xtype = repmat('C',length(lb),1);
xtype(model.binary_variables) = 'B';
xtype(model.integer_variables) = 'I';
x0 = model.x0;
opts = model.options.scip;
switch model.options.verbose
    case 0
        opts.display = 'off';
    otherwise
        opts.display = 'iter';
end
if model.options.savedebug    
    save scipnldebug obj con A ru rl cl cu xtype lb ub x0 opts
end

solvertime = tic;
[x,fval,exitflag,info] = opti_scipnl(obj,A,rl,ru,lb,ub,con,cl,cu,xtype,[],opts);%,x0,opts);
solvertime = toc(solvertime);

% Check, currently not exhaustive...
switch exitflag
    case 0
        if ~isempty(strfind(info.Status,'Exceed'))
            problem = 3;
        else
            problem = 9;
        end
    case 1
        problem = 0;
    case {2,-1}
        problem = 1;
    case 3
        problem = 2;
    case {4,5}
        problem = 11;
    otherwise
        problem = 9;
end

% Save all data sent to solver?
if model.options.savesolverinput
    solverinput.obj = obj;
    solverinput.A = A;
    solverinput.rl = rl;
    solverinput.ru = ru;
    solverinput.lb = lb;
    solverinput.ub = ub;
    solverinput.con = con;
    solverinput.cl = cl;
    solverinput.cu = cu;
    solverinput.xtype = xtype;
    solverinput.opts = opts;
else
    solverinput = [];
end

% Save all data from the solver?
if model.options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fval = fval;
    solveroutput.exitflag = exitflag;
    solveroutput.info=info;
    solveroutput.allsol=allsol;
else
    solveroutput = [];
end

if ~isempty(x)
    x = RecoverNonlinearSolverSolution(model,x);
end

% Standard interface
output = createoutput(x,[],[],problem,'SCIP-NL',solverinput,solveroutput,solvertime);



function z = createOneterm(i,model)

[evalPos,pos] = ismember(i,model.evalVariables);
if pos
    % This is a nonlinear operator
    map = model.evalMap{pos};
    % We might hqave vectorized things to compute several expressions at the same time
    if length(map.computes) == 1 && length(map.variableIndex) > 1
        if isequal(map.fcn,'plog')
            j1 = find(model.linearindicies == map.variableIndex(1));
            j2 = find(model.linearindicies == map.variableIndex(2));
            if isempty(j1)
                z1 = createmonomstring(model.monomtable(map.variableIndex(1),:),model);
            else
                z1 = ['x(' num2str(j1) ')'];
            end
            if isempty(j2)
                z1 =  createmonomstring(model.monomtable(map.variableIndex(2),:),model);
            else
                z2 =  ['x(' num2str(j2) ')'];
            end
            z = [z1 '*log(' z2 '/' z1 ')'];
        else
        z = [map.fcn '('];
        for j = map.variableIndex
            jl = find(model.linearindicies == j);
            if isempty(jl)
                z =  [z createmonomstring(model.monomtable(j,:),model)];
            else
                z =  [z 'x(' num2str(jl) ')'];
            end
            if j == length(map.variableIndex)
                z = [z ')'];
            else
                z = [z ','];
            end
        end
        end
    else
        j = find(map.computes == i);
        j = map.variableIndex(j);
        % we have f(x(j)), but have to map back to linear indicies
        jl = find(model.linearindicies == j);
        if isempty(jl)
            z =  [map.fcn '(' createmonomstring(model.monomtable(j,:),model)  ')'];
        else
            z =  [map.fcn '(x(' num2str(jl) '))'];
        end
    end
else
    i = find(model.linearindicies == i);
    z =  ['x(' num2str(i) ')'];
end

function monomstring = createmonomstring(v,model)

monomstring = '';
for i = find(v)
    z = createOneterm(i,model);
    if v(i)==1
        monomstring = [monomstring z '*'];
    else
        monomstring = [monomstring z '^' num2str(v(i),12) '*'];
    end
end
monomstring = monomstring(1:end-1);


function string = createmodelstring(row,model)
string = '';
index = find(row);
for i = index(:)'
    monomstring = createmonomstring(model.monomtable(i,:),model);
    string = [string  num2str(row(i),12) '*' monomstring '+'];
end
string = string(1:end-1);
string = strrep(string,'+-','-');
string = strrep(string,')-1*',')-');
string = strrep(string,')+1*',')+');

function string = createQstring(Q,model)
% Special code for the quadratic term
string = '';
[ii,jj,c]= find(triu(Q));
for i = 1:length(ii)
    if ii(i)==jj(i)
        % Quadratic term
        monomstring = createmonomstring(sparse(1,ii(i),1,1,size(Q,1)),model);       
        string = [string  num2str(c(i),12) '*' monomstring '^2' '+'];
    else
        % Bilinear term
        monomstring1 = createmonomstring(sparse(1,ii(i),1,1,size(Q,1)),model);
        monomstring2 = createmonomstring(sparse(1,jj(i),1,1,size(Q,1)),model);
        string = [string  num2str(c(i),12) '*2*' monomstring1 '*' monomstring2 '+'];
    end
end
string = string(1:end-1);
string = strrep(string,'+-','-');
string = strrep(string,'-1*','-');
string = strrep(string,'+1*','+');


