function diagnostic = solvesdp_multiple(varargin);

yalmiptime = clock;
h = varargin{2};
varargin{2} = sum(recover(depends(h)));

if ~is(h,'linear')
    error('Parts of your matrix objective is not linear (multiple solutions can currently only be obtained for linear objectives)');
end

if is(h,'complex')
    error('Parts of your matrix objective is complex-valued (which makes no sense since complex numbers have no natural ordering');
end

if nargin<3
    ops = sdpsettings('saveyalmipmodel',1);
    varargin{3} = ops;
else
    ops = varargin{3};
    varargin{3}.saveyalmipmodel = 1;
end
varargin{3}.pureexport = 1;
yalmiptime = clock;
model = solvesdp(varargin{:});
diagnostic.yalmiptime = etime(clock,yalmiptime);
yalmiptime = clock;

h_variables = getvariables(h);
h_base = getbase(h);
for i = 1:length(h)
    model.f = h_base(i,1);
    model.c = mapObjective(h_variables,model.used_variables,h_base(i,2:end));
    
    % *************************************************************************
    % TRY TO SOLVE PROBLEM
    % *************************************************************************
    if ops.debug
        eval(['output = ' model.solver.call '(model);']);
    else
        try
            eval(['output = ' model.solver.call '(model);']);
        catch
            output.Primal = zeros(length(model.c),1)+NaN;
            output.Dual  = [];
            output.Slack = [];
            output.solvertime   = nan;
            output.solverinput  = [];
            output.solveroutput = [];
            output.problem = 9;
            output.infostr = yalmiperror(output.problem,lasterr);
        end
    end
    
    if ops.dimacs
        try
            b = -model.c;
            c = model.F_struc(:,1);
            A = -model.F_struc(:,2:end)';
            x = output.Dual;
            y = output.Primal;
            % FIX this nonlinear crap (return variable type in
            % compilemodel)
            if options.relax == 0 & any(full(sum(model.monomtable,2)~=0))
                if ~isempty(find(sum(model.monomtable | model.monomtable,2)>1))
                    z=real(exp(model.monomtable*log(y+eps)));
                    y = z;
                end
            end
            
            if isfield(output,'Slack')
                s = output.Slack;
            else
                s = [];
            end
            
            dimacs = computedimacs(b,c,A,x,y,s,model.K);
        catch
            dimacs = [nan nan nan nan nan nan];
        end
    else
        dimacs = [nan nan nan nan nan nan];
    end
    
    % ********************************
    % ORIGINAL COORDINATES
    % ********************************    
    if isempty(output.Primal)
       output.Primal = zeros(size(model.recoverdata.H,2),1);
    end
    output.Primal(:,i) = model.recoverdata.x_equ+model.recoverdata.H*output.Primal;
    
    % ********************************
    % OUTPUT
    % ********************************
    diagnostic.yalmiptime = diagnostic.yalmiptime + etime(clock,yalmiptime)-output.solvertime;
    diagnostic.solvertime(:,i) = output.solvertime;    
    try
        diagnostic.info = output.infostr;
    catch
        diagnostic.info = yalmiperror(output.problem,model.solver.tag);
    end
    diagnostic.problem(:,i) = output.problem;
    diagnostic.dimacs(:,i) = dimacs;
    
    % Some more info is saved internally
    solution_internal = diagnostic;
    solution_internal.variables = model.recoverdata.used_variables(:);
    solution_internal.optvar = output.Primal(:,i);
    
    if ~isempty(model.parametric_variables)
        diagnostic.mpsol = output.solveroutput;
        options.savesolveroutput = actually_save_output;
    end;
    
    if model.options.savesolveroutput
        diagnostic.solveroutput = output.solveroutput;
    end
    if model.options.savesolverinput
        diagnostic.solverinput = output.solverinput;
    end
    if model.options.saveyalmipmodel
        diagnostic.yalmipmodel = model;
    end
    
    if ops.warning && warningon && isempty(findstr(output.infostr,'No problems detected'))
        disp(['Warning: ' output.infostr]);
    end
    
    if ismember(output.problem,ops.beeponproblem)
        try
            beep; % does not exist on all ML versions
        catch
        end
    end
    
    % And we are done! Save the result
    if ~isempty(output.Primal)
        yalmip('setsolution',solution_internal,i);
    end
    if model.options.saveduals & model.solver.dual
        if isempty(model.Fremoved) | (nnz(model.Q)>0)
            try
                setduals(F,output.Dual,model.K);
            catch
            end
        else
            try
                % Duals related to equality constraints/free variables
                % have to be recovered b-A*x-Ht == 0
                b = -model.oldc;
                A = -model.oldF_struc(1+model.oldK.f:end,2:end)';
                H = -model.oldF_struc(1:model.oldK.f,2:end)';
                x = output.Dual;
                b_equ = b-A*x;
                newdual =  H\b_equ;
                setduals(model.Fremoved + F,[newdual;output.Dual],model.oldK);
            catch
                % this is a new feature...
                disp('Dual recovery failed. Please report this issue.');
            end
        end
    end       
end


function newbase = mapObjective(local_vars,global_vars,base)
newbase = spalloc(length(global_vars),1,nnz(base));
for i = 1:length(local_vars)
    j = find(local_vars(i)==global_vars);
    newbase(j) = base(i);
end

function yesno = warningon

s = warning;
yesno = isequal(s,'on');


