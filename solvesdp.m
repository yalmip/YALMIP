function diagnostic = solvesdp(varargin)
%SOLVESDP Computes solution to optimization problem
%
%   DIAGNOSTIC = SOLVESDP(F,h,options) is the common command to
%   solve optimization problems of the following kind
%
%    min        h
%    subject to
%            F >=(<=,==) 0
%
%   NOTES
%    Despite the name, SOLVESDP is the interface for solving all
%    supported problem classes (LP, QP, SOCP, SDP, BMI, MILP, MIQP,...)
%
%    To obtain solution for a variable, use DOUBLE.
%
%    To obtain dual variable for a constraint, use DUAL.
%
%    See YALMIPERROR for error codes returned in output.
%
%   OUTPUT
%     diagnostic : Diagnostic information
%
%   INPUT
%     F          : Object describing the constraints. Can be [].
%     h          : SDPVAR object describing the objective h(x). Can be [].
%     options    : Options structure. See SDPSETTINGS. Can be [].
%
%   EXAMPLE
%    A = randn(15,5);b = rand(15,1)*5;c = randn(5,1);
%    x = sdpvar(5,1);
%    solvesdp(set(A*x<=b),c'*x);double(x)
%
%   See also DUAL, @SDPVAR/DOUBLE, SDPSETTINGS, YALMIPERROR

% Author Johan Löfberg
% $Id: solvesdp.m,v 1.75 2010-04-06 06:32:49 joloef Exp $

yalmiptime = clock; % Let us see how much time we spend

% Avoid warning
if length(varargin)>=2
    if isa(varargin{2},'double')
        varargin{2} = [];
    end
end

if length(varargin)>=2
    if isa(varargin{2},'sdpvar') && prod(size(varargin{2}))>1
        % Several objectives
        diagnostic = solvesdp_multiple(varargin{:});
        return
    end
end

% Arrrgh, new format with logdet much better, but we have to
% take care of old code, requires some testing...
varargin = combatible({varargin{:}});
nargin = length(varargin);

% *********************************
% CHECK INPUT
% *********************************


if nargin<1
    help solvesdp
    return
else
    F = varargin{1};
    if isa(F,'constraint')
        F = set(F);
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
            Fnew = Fnew + set(true(F(i)));
        end
        if isnan(Fnew)
            error('First argument (F) should be a constraint object.');
        else
            F = Fnew;
        end
    elseif isempty(F)
        F = set([]);
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
        error('Second argument (the objective function h) should be an sdpvar or logdet object (or empty).');
    end
    if isa(h,'logdet')
        logdetStruct.P     = getP(h);
        logdetStruct.gain  = getgain(h);
%         if any(logdetStruct.gain>0)
%             warning('Nonconvex terms! Perhaps you mean -logdet(P)...')
%             diagnostic.yalmiptime = etime(clock,yalmiptime);
%             diagnostic.solvertime = 0;
%             diagnostic.info = yalmiperror(-2,'YALMIP');
%             diagnostic.problem = -2;
%             return
%         end
        h = getcx(h);
        if isempty(F)
            F = set([]);
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
if length(F) > 0
    rand_declarations = is(F,'random');
    if any(rand_declarations)
    %    diagnostic = solverandom(F(find(~rand_declarations)),h,options,recover(getvariables(sdpvar(F(find(unc_declarations))))));
        return
    end
end


% Call robust solver?
if length(F) > 0
    unc_declarations = is(F,'uncertain');
    if any(unc_declarations)
        diagnostic = solverobust(F(find(~unc_declarations)),h,options,recover(getvariables(sdpvar(F(find(unc_declarations))))));
        return
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
            %im = im-diag(diag(im));
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

% ******************************************
% DID WE SELECT THE BMILIN SOLVER (obsolete)
% ******************************************
if strcmpi(solver.tag,'bmilin')
    diagnostic = callbmilin(F,h,options);
    return
end

% ******************************************
% DID WE SELECT THE BMIALT SOLVER (obsolete)
% ******************************************
if strcmp(solver.tag,'bmialt')
    diagnostic = callbmialt(F,h,options);
    return
end

%******************************************
% DID WE SELECT THE MPT solver (backwards comb)
%******************************************
actually_save_output = interfacedata.options.savesolveroutput;
if strcmpi(solver.tag,'mpt-2') | strcmpi(solver.tag,'mpt-3') | strcmpi(solver.tag,'mpcvx') | strcmpi(solver.tag,'mplcp')    
    interfacedata.options.savesolveroutput = 1;
    if isempty(interfacedata.parametric_variables)
        if (nargin < 4 | ~isa(varargin{4},'sdpvar'))
            error('You must specify parametric variables.')
        else
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

% *************************************************************************
% TRY TO SOLVE PROBLEM
% *************************************************************************
if options.debug
    eval(['output = ' solver.call '(interfacedata);']);
else
    try
        eval(['output = ' solver.call '(interfacedata);']);
    catch
        output.Primal = zeros(length(interfacedata.c),1)+NaN;
        output.Dual  = [];
        output.Slack = [];
        output.solvertime   = nan;
        output.solverinput  = [];
        output.solveroutput = [];
        output.problem = 9;
        output.infostr = yalmiperror(output.problem,lasterr);
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

if options.warning & warningon & isempty(findstr(diagnostic.info,'No problems detected'))
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
                Ktemp.l = length(output.qcDual);
                Ktemp.f = 0;
                Ktemp.q = 0;
                Ktemp.s = 0;
                Ktemp.r = 0;
                setduals(ForiginalQuadratics,-output.qcDual,Ktemp);
            end
        end
    end
end

function newinputformat = combatible(varargin)

varargin = varargin{1};

classification = 0;
% 0 : Ambigious
% 1 : Old
% 2 : New

% Try some fast methods to determine...
m = length(varargin);
if m==1
    classification = 2;
elseif m>=3 && isstruct(varargin{3})
    classification = 2;
elseif m>=4 && isstruct(varargin{4})
    classification = 1;
elseif m>=2 && isa(varargin{2},'lmi')
    classification = 1;
elseif m>=3 && isa(varargin{3},'sdpvar')
    classification = 1;
elseif m>=2 && isa(varargin{2},'sdpvar') & min(size(varargin{2}))==1
    classification = 2;
elseif m>=2 && isa(varargin{2},'sdpvar') & prod(size(varargin{2}))>=1
    classification = 1;
elseif m>=2 && isa(varargin{2},'logdet')
    classification = 2;
elseif m==2 && isempty(varargin{2})
    classification = 2;
end

if classification==0
    warning('I might have interpreted this problem wrong due to the new input format in version 3. To get rid of this warning, use an options structure');
    classification = 2;
end

if classification==2
    newinputformat = varargin;
else
    newinputformat = varargin;
    P = varargin{2};
    % 99.9% of the cases....
    if isempty(P)
        newinputformat = {newinputformat{[1 3:end]}};
    else
        if isa(P,'lmi')
            P = sdpvar(P);
        end
        if m>=3
            cxP = newinputformat{3}-logdet(P);
            newinputformat{3}=cxP;
        else
            cxP = -logdet(P);
            newinputformat{3}=cxP;
        end
        newinputformat = {newinputformat{[1 3:end]}};
    end
end

function yesno = warningon

s = warning;
yesno = isequal(s,'on');