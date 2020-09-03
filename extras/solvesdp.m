function diagnostic = solvesdp(varargin)
%SOLVESDP Obsolete command, please use OPTIMIZE

yalmiptime = clock; % Let us see how much time we spend

% *********************************
% CHECK INPUT
% *********************************
nargin = length(varargin);

% First check of objective for early transfer to multiple solves
if nargin>=2
    if isa(varargin{2},'double')
        varargin{2} = [];
    elseif isa(varargin{2},'sdpvar') && numel(varargin{2})>1
        % Several objectives
        diagnostic = solvesdp_multiple(varargin{:});
        return
    end
end

% Early jump to bisection solver
if nargin >= 3
    if isa(varargin{3},'struct')
        if isfield(varargin{3},'solver')
            if isequal(varargin{3}.solver,'bisection')    
                if any(is(varargin{1},'sos'))
                    [F_sos,h_sos] = compilesos(varargin{1},varargin{2},varargin{3});
                    varargin{1} = F_sos;
                    varargin{2} = h_sos;
                end
                diagnostic = bisection(varargin{:});
                return
            end
        end
    end
end

if nargin<1
    help solvesdp
    return
else
    F = varargin{1};
    if isa(F,'constraint')
        F = lmi(F);
    end
    if isa(F,'lmi')
        F = flatten(F);
    end
    
    if isa(F,'sdpvar')
        % We do allow sloppy coding of logic constraints, i.e writing a
        % constraints as [a|b true(a)]
        Fnew = [];
        for i = 1:length(F)
            if length(getvariables(F(i)))>1
                Fnew = nan;
                break
            end
            operator = yalmip('extstruct',getvariables(F(i)));
            if isempty(operator)
                Fnew = nan;
                break
            end
            if length(operator)>1
                Fnew = nan;
                break
            end
            if ~strcmp(operator.fcn,'or')
                Fnew = nan;
                break
            end
            Fnew = Fnew + (true(F(i)));
        end
        if isnan(Fnew)
            error('First argument (F) should be a constraint object.');
        else
            F = Fnew;
        end
    elseif isempty(F)
        F = lmi([]);
    elseif ~isa(F,'lmi')
        error('First argument (F) should be a constraint object.');      
    end
end

if nargin>=2
    h = varargin{2};
    if isa(h,'double')
        h = [];
    end
    if ~(isempty(h) | isa(h,'sdpvar') | isa(h,'logdet') |  isa(h,'ncvar'))
        if isa(h,'struct')
            error('Second argument (the objective function h) should be an sdpvar or logdet object (or empty). It appears as if you sent an options structure in the second argument.');
        else
            error('Second argument (the objective function h) should be an sdpvar or logdet object (or empty).');
        end
    end
    if isa(h,'logdet')
        logdetStruct.P     = getP(h);
        logdetStruct.gain  = getgain(h);
        h = getcx(h);
        if isempty(F)
            F = ([]);
        end
    else
        logdetStruct = [];
    end
else
    logdetStruct = [];
    h = [];   
end

if ~isempty(F)
    if any(is(F,'sos'))        
        diagnostic = solvesos(varargin{:});
        return
    end
end

if isa(h,'sdpvar')
    if is(h,'complex')
        error('Complex valued objective does not make sense.');
    end
end
    
if nargin>=3
    options = varargin{3};
    if ~(isempty(options) | isa(options,'struct'))
        error('Third argument (options) should be an sdpsettings struct (or empty).');
    end
    if isempty(options)
        options = sdpsettings;
    end
else
    options = sdpsettings;
end
options.solver = lower(options.solver);

% If user has logdet term, but no preference on solver, we try to hook up
% with SDPT3 if possible.
if ~isempty(logdetStruct)
    if strcmp(options.solver,'')
     %   options.solver = 'sdpt3,*';
    end
end

% Call chance solver?
% if length(F) > 0
%     rand_declarations = is(F,'random');
%     if any(rand_declarations)
%     %    diagnostic = solverandom(F(find(~rand_declarations)),h,options,recover(getvariables(sdpvar(F(find(unc_declarations))))));
%         return
%     end
% end


% Call robust solver?
if length(F) > 0
    unc_declarations = is(F,'uncertain');
    if any(unc_declarations)
        try
            diagnostic = solverobust(F(find(~unc_declarations)),h,options,recover(getvariables(sdpvar(F(find(unc_declarations))))));
            return
        catch
            if strfind(lasterr,'Undefined function ''solverobust'' ') 
                display('***');
                display('It looks like you have failed to add yalmip/modules/robust to your path');
                display('Information on required paths https://yalmip.github.io/tutorial/installation/');
                display(['For complete path, use addpath(genpath(''' fileparts(which('yalmiptest.m')) '''))']);
                display('***');
                error(lasterr)
            else
                error(lasterr)
            end
        end        
    end
end
    
if isequal(options.solver,'mpt') | nargin>=4
    solving_parametric = 1;
else
    solving_parametric = 0;
end
    
% Just for safety
if isempty(F) & isempty(logdetStruct)
    F = lmi;
end

if any(is(F,'sos'))
    error('You have SOS constraints. Perhaps you meant to call SOLVESOS.');
end

% Super stupido
if length(F) == 0 & isempty(h) & isempty(logdetStruct)
   diagnostic.yalmiptime = 0;
   diagnostic.solvertime = 0;
   diagnostic.info = 'No problems detected (YALMIP)';
   diagnostic.problem = 0;
   diagnostic.dimacs = [NaN NaN NaN NaN NaN NaN];
   return
end

% Dualize the problem?
if ~isempty(F)
    if options.dualize == -1
        sdp = find(is(F,'sdp'));
        if ~isempty(sdp)
            if all(is(F(sdp),'sdpcone'))
                options.dualize = 1;
            end
        end
    end
end
if options.dualize == 1   
    [Fd,objd,aux1,aux2,aux3,complexInfo] = dualize(F,h,[],[],[],options);
    options.dualize = 0;
    diagnostic = solvesdp(Fd,-objd,options);
    if ~isempty(complexInfo)
        for i = 1:length(complexInfo.replaced)
            n = size(complexInfo.replaced{i},1);
            re = 2*double(complexInfo.new{i}(1:n,1:n));            
            im = 2*double(complexInfo.new{i}(1:n,n+1:end));
            im=triu((im-im')/2)-(triu((im-im')/2))';
            assign(complexInfo.replaced{i},re + sqrt(-1)*im);
        end
    end
    return
end

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DID WE SELECT THE MOMENT SOLVER
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
if isequal(options.solver,'moment')
    if ~isempty(logdetStruct)
        error('Cannot dualize problems with logaritmic objective')
    end
    options.solver = options.moment.solver;
    [diagnostic,x,momentdata] = solvemoment(F,h,options,options.moment.order);
    diagnostic.momentdata = momentdata;
    diagnostic.xoptimal = x;
    return
end

% ******************************************
% COMPILE IN GENERALIZED YALMIP FORMAT
% ******************************************
[interfacedata,recoverdata,solver,diagnostic,F,Fremoved,ForiginalQuadratics] = compileinterfacedata(F,[],logdetStruct,h,options,0,solving_parametric);

% ******************************************
% FAILURE?
% ******************************************
if ~isempty(diagnostic)
    diagnostic.yalmiptime = etime(clock,yalmiptime);
    return
end

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DID WE SELECT THE LMILAB SOLVER WITH A KYP
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
if  strcmpi(solver.tag,'lmilab') & any(is(F,'kyp'))
    [diagnostic,failed] = calllmilabstructure(F,h,options);
    if ~failed % Did this problem pass (otherwise solve using unstructured call)
        diagnostic.yalmiptime = etime(clock,yalmiptime)-diagnostic.solvertime;
        return
    end
end

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DID WE SELECT THE KYPD SOLVER
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
if strcmpi(solver.tag,'kypd')
    diagnostic = callkypd(F,h,options);
    diagnostic.yalmiptime = etime(clock,yalmiptime)-diagnostic.solvertime;
    return
end

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DID WE SELECT THE STRUL SOLVER
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
if strfind(solver.tag,'STRUL')
    diagnostic = callstrul(F,h,options);
    diagnostic.yalmiptime = etime(clock,yalmiptime)-diagnostic.solvertime;
    return
end

%******************************************
% DID WE SELECT THE MPT solver (backwards comb)
%******************************************
actually_save_output = interfacedata.options.savesolveroutput;
if strcmpi(solver.tag,'mpt-2') | strcmpi(solver.tag,'mpt-3') | strcmpi(solver.tag,'mpcvx') | strcmpi(solver.tag,'mplcp') | strcmpi(solver.tag,'pop')    
    interfacedata.options.savesolveroutput = 1;
    if isempty(interfacedata.parametric_variables)
        if (nargin < 4 | ~isa(varargin{4},'sdpvar'))
            error('You must specify parametric variables.')
        else
            varargin{4} = reshape(varargin{4},[],1);
            interfacedata.parametric_variables = [];
            for i = 1:length(varargin{4})
                  interfacedata.parametric_variables = [interfacedata.parametric_variables;find(ismember(recoverdata.used_variables,getvariables(varargin{4}(i))))];
            end            
            if isempty(varargin{5})
                interfacedata.requested_variables = [];
            else
                interfacedata.requested_variables = [];
                for i = 1:length(varargin{5})
                    interfacedata.requested_variables = [interfacedata.requested_variables;find(ismember(recoverdata.used_variables,getvariables(varargin{5}(i))))];
                end
            end
        end
    end
end

% *************************************************************************
% Just return the YALMIP model. Used when solving multiple objectives
% *************************************************************************
if isfield(options,'pureexport')
    interfacedata.recoverdata = recoverdata;
    diagnostic = interfacedata;        
    return
end

if strcmpi(solver.version,'geometric') || (strcmpi(solver.tag,'bnb') && strcmpi(solver.lower.version,'geometric'))
    % Actual linear user variables
    if options.assertgpnonnegativity
        check = find(interfacedata.variabletype==0);
        check = setdiff(check,interfacedata.aux_variables);
        check = setdiff(check,interfacedata.evalVariables);
        check = setdiff(check,interfacedata.extended_variables);
        [lb,ub] = findulb(interfacedata.F_struc,interfacedata.K);
        if ~all(lb(check)>=0)
            % User appears to have explictly selected a GP solver
            userdirect = ~isempty(strfind(options.solver,'geometric')) || ~isempty(strfind(options.solver,'mosek')) || ~isempty(strfind(options.solver,'gpposy'));
            userindirect = strcmpi(solver.tag,'bnb') && strcmpi(solver.lower.version,'geometric');
            if userdirect || userindirect
                % There are missing non-negativity bounds
                output = createOutputStructure(zeros(length(interfacedata.c),1)+NaN,[],[],18,yalmiperror(18,''),[],[],nan);
                diagnostic.yalmiptime = etime(clock,yalmiptime);
                diagnostic.solvertime = output.solvertime;
                try
                    diagnostic.info = output.infostr;
                catch
                    diagnostic.info = yalmiperror(output.problem,solver.tag);
                end
                diagnostic.problem = output.problem;
                if options.dimacs
                    diagnostic.dimacs = dimacs;
                end
                return
            else
                % YALMIP selected solver and picked a GP solver. As this is
                % no GP, we call again, but this time explicitly tell
                % YALMIP that it isn't a GP
                options.thisisnotagp = 1;
                varargin{3} = options;
                diagnostic = solvesdp(varargin{:});
                return
            end
        end
    end
end

% *************************************************************************
% TRY TO SOLVE PROBLEM
% *************************************************************************
if options.debug
    eval(['output = ' solver.call '(interfacedata);']);
else
    try
        eval(['output = ' solver.call '(interfacedata);']);
    catch
        output = createOutputStructure(zeros(length(interfacedata.c),1)+NaN,[],[],9,yalmiperror(9,lasterr),[],[],nan);        
    end
end

if options.dimacs
    try       
        b = -interfacedata.c;
        c = interfacedata.F_struc(:,1);
        A = -interfacedata.F_struc(:,2:end)';
        x = output.Dual;
        y = output.Primal;
        % FIX this nonlinear crap (return variable type in
        % compileinterfacedata)
        if options.relax == 0 & any(full(sum(interfacedata.monomtable,2)~=0))
            if ~isempty(find(sum(interfacedata.monomtable | interfacedata.monomtable,2)>1))                
                z=real(exp(interfacedata.monomtable*log(y+eps)));                
                y = z;
            end
        end
        
        if isfield(output,'Slack')
            s = output.Slack;
        else
            s = [];
        end
                    
        dimacs = computedimacs(b,c,A,x,y,s,interfacedata.K);
    catch
        dimacs = [nan nan nan nan nan nan];
    end
else
    dimacs = [nan nan nan nan nan nan];
end

% ********************************
% ORIGINAL COORDINATES
% ********************************
output.Primal = recoverdata.x_equ+recoverdata.H*output.Primal;

% ********************************
% OUTPUT
% ********************************
diagnostic.yalmipversion = yalmip('ver');
diagnostic.matlabversion = version;
diagnostic.yalmiptime = etime(clock,yalmiptime)-output.solvertime;
diagnostic.solvertime = output.solvertime;
try
    diagnostic.info = output.infostr;
catch   
    diagnostic.info = yalmiperror(output.problem,solver.tag);
end
diagnostic.problem = output.problem;
if options.dimacs
    diagnostic.dimacs = dimacs;
end

% Some more info is saved internally
solution_internal = diagnostic;
solution_internal.variables = recoverdata.used_variables(:);
solution_internal.optvar = output.Primal;

if ~isempty(interfacedata.parametric_variables)
    diagnostic.mpsol = output.solveroutput;
    options.savesolveroutput = actually_save_output;
end;

if interfacedata.options.savesolveroutput
    diagnostic.solveroutput = output.solveroutput;
end
if interfacedata.options.savesolverinput
    diagnostic.solverinput = output.solverinput;
end
if interfacedata.options.saveyalmipmodel
    diagnostic.yalmipmodel = interfacedata;
end

if options.warning & warningon & isempty(strfind(diagnostic.info,'No problems detected'))
    disp(['Warning: ' output.infostr]);
end

if ismember(output.problem,options.beeponproblem)
    try
        beep; % does not exist on all ML versions
    catch
    end
end

% And we are done! Save the result
if ~isempty(output.Primal)
    if size(output.Primal,2)>1
        for j = 1:size(output.Primal,2)
            temp = solution_internal;
            temp.optvar = temp.optvar(:,j);
            yalmip('setsolution',temp,j);
        end
    else
        yalmip('setsolution',solution_internal);
    end
end
if interfacedata.options.saveduals & solver.dual 
    if isempty(interfacedata.Fremoved) | (nnz(interfacedata.Q)>0)
        try
            setduals(F,output.Dual,interfacedata.K);
        catch
        end
    else
        try
            % Duals related to equality constraints/free variables
            % have to be recovered b-A*x-Ht == 0
            b = -interfacedata.oldc;          
            A = -interfacedata.oldF_struc(1+interfacedata.oldK.f:end,2:end)';
            H = -interfacedata.oldF_struc(1:interfacedata.oldK.f,2:end)';
            x = output.Dual;
            b_equ = b-A*x;
            newdual =  H\b_equ;
            setduals(interfacedata.Fremoved + F,[newdual;output.Dual],interfacedata.oldK);
        catch
             % this is a new feature...
            disp('Dual recovery failed. Please report this issue.');
        end
    end
end
% Hack to recover original QCQP duals from gurobi
if strcmp(solver.tag,'GUROBI-GUROBI')
    if length(ForiginalQuadratics) > 0
        if isfield(output,'qcDual')
            if length(output.qcDual) == length(ForiginalQuadratics)               
             %   Ktemp.l = length(output.qcDual);
             %   Ktemp.f = 0;
             %   Ktemp.q = 0;
             %   Ktemp.s = 0;
             %   Ktemp.r = 0;
                Ftemp = F + ForiginalQuadratics;
                K = interfacedata.K;
                Ktemp = K;
                Ktemp.l = Ktemp.l + length(ForiginalQuadratics);
                tempdual = output.Dual;
                tempdual = [tempdual(1:K.f + K.l);-output.qcDual;tempdual(1+K.f+K.l:end)];
                setduals(Ftemp,tempdual,Ktemp);
%                setduals(ForiginalQuadratics,-output.qcDual,Ktemp);
            end
        end
    end
end

function yesno = warningon

s = warning;
yesno = isequal(s,'on');
