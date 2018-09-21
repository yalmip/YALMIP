function [model,recoverdata,diagnostic,interfacedata] = export(varargin)
%EXPORT  Exports YALMIP problem to solver specific format
%
%   [MODEL,RECOVERYMODEL,DIAGNOSTIC,INTERNAL] = EXPORT(F,h,options) is the
%   common command to export optimization problems of the following kind
%
%    min        h
%    subject to
%            F >(=) 0
%
%
%   The MODEL is exported in the format defined by the solver chosen
%   in the options structure, or automatically chosen by YALMIP.
%
%   If the solver format not is support by EXPORT,the YALMIP model used to
%   call the solver is returned)
%
%   If YALMIP by some reason failed to generate a model, the DIAGNOSTIC
%   variable will be non-empty.
%
%   The fourth output is the internal model used by YALMIP to communicate
%   with the generalized solver interfaces.
%
%   The RECOVERYMODEL is used to relate a solution of the exported model
%   to the original variables in YALMIP.

nargin = length(varargin);
% *********************************
% CHECK INPUT
% *********************************
if nargin<1
    help export
    return
else
    F = varargin{1};
    % Check for wrong syntax
    if ~isempty(F) & ~isa(F,'lmi') & ~isa(F,'constraint')
        error('First argument should be a SET object')
    end

    if isa(F,'constraint')
        F = lmi(F);
    end
end

model = [];
recoverdata = [];
diagnostic = [];
interfacedata = [];

if nargin>=2
    h = varargin{2};
    if isa(h,'double')
        h = [];
    end
    if ~(isempty(h) | isa(h,'sdpvar') | isa(h,'logdet'))
        error('Second argument (the objective function h) should be an sdpvar or logdet object (or empty).');
    end
    if isa(h,'logdet')
        logdetStruct.P  = getP(h);
        logdetStruct.gain  = getgain(h);
        if any(logdetStruct.gain>0)
            warning('Perhaps you mean -logdet(P)...')
            diagnostic.yalmiptime = etime(clock,yalmiptime);
            diagnostic.solvertime = 0;
            diagnostic.info = yalmiperror(-2,'YALMIP');
            diagnostic.problem = -2;
            return
        end
        h = getcx(h);
        if isempty(F)
           F = ([]);
        end

    else
      logdetStruct = [];
    end
else
    h = [];
    logdetStruct = [];
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

if nargin<6
    if isequal(options.solver,'')
        findallsolvers = 1;
    else
        findallsolvers = 0;
    end
else
    findallsolvers = varargin{6};
end

% Just for safety
if isempty(F) & isempty(logdetStruct)
    F = lmi;
end

if any(is(F,'uncertain'))
    [F,h] = robustify(F,h,options);
end

F = flatten(F);

% ******************************************
% Export SOS problem to SOS first
% ******************************************
 if any(is(F,'sos'))
     old =  options.verbose;
     options.verbose = max(options.verbose - 1,0);
     [F,h] = compilesos(F,h,options);
     options.verbose = old;
 end

% ******************************************
% COMPILE IN GENERALIZED YALMIP FORMAT
% ******************************************
if ~isempty(F) & any(is(F,'parametric'))
    % Special code, do not mix with standard case, in case there is some
    % bug (this has not been tested)
    [interfacedata,recoverdata,solver,diagnostic,F] = compileinterfacedata(F,[],logdetStruct,h,options,findallsolvers,1);
else
    [interfacedata,recoverdata,solver,diagnostic,F] = compileinterfacedata(F,[],logdetStruct,h,options,findallsolvers);
end

if ~isempty(diagnostic)
    model = [];
    recoverdata = [];
    return
end

% Not official yet
if nargin == 5
    model=interfacedata;
    return
end

% ******************************************
% CONVERT
% ******************************************
switch lower(solver.tag)

    case 'cplex-ibm'

        % Hack to handle CPLEX slow treatment of parameters. Remove all
        % default settings, so optimizer runs fast
        interfacedata.options = prunecplexoptions(interfacedata.options);
        model = yalmip2cplex(interfacedata);

    case 'ecos'
        model = yalmip2ecos(interfacedata);

    case 'osqp'
        model = yalmip2osqp(interfacedata);

    case 'cbc'
        model = yalmip2cbc(interfacedata);

    case 'dsdp-opti'
        model = yalmip2optidsdp(interfacedata);

    case 'gurobi-gurobi'
        model = yalmip2gurobi(interfacedata);

    case 'gurobi-mex'
        model = yalmip2gurobimex(interfacedata);

    case 'cplex-cplexint'
        [model.H,model.C,model.A,model.B,model.LB,model.UB,model.QC,model.VARTYPE,model.INDEQ,model.PARAM,model.OPTIONS] = cplex2yalmip(interfacedata);

    case {'mosek-socp','mosek-lp/qp','mosek-geometric','mosek-sdp'}
        if interfacedata.K.s(1)>0
            model.prob = yalmip2SDPmosek(interfacedata);
        else
            model.prob = yalmip2mosek(interfacedata);
        end
        model.param = interfacedata.options.mosek;

    case 'linprog'
        model = yalmip2quadprog(interfacedata);
        model = rmfield(model,'Q');

    case 'quadprog'
        model = yalmip2quadprog(interfacedata);

    case 'intlinprog'
        model = yalmip2intlinprog(interfacedata);

    case {'sedumi-1.05','sedumi-1.1','sedumi-1.3'}
        model = yalmip2sedumi(interfacedata);

   case {'scs-direct','scs-indirect'}
        model = yalmip2scs(interfacedata);

    case {'powersolver'}
        model = yalmip2powersolver(interfacedata);

    case 'csdp'
        model = yalmip2csdp(interfacedata);

    case 'dsdp-5'
        model = yalmip2dsdp(interfacedata);

    case 'sdpa-m'
        model = yalmip2sdpa(interfacedata);

    case {'sdpt3-3.1','sdpt3-4'}
        % Convert from internal (sedumi-like) format
        if isequal(interfacedata.K.m,0)
           model = yalmip2sdpt3(interfacedata);
        else
            error('MAXDET models still not supported in SDPT3 export')
        end

    case {'glpk-glpkmex','glpk-glpkmex-cc'}
        model = yalmip2glpkmex(interfacedata);

    case 'pensdp-penopt'
        model = yalmip2pensdp(interfacedata);

    case 'penlab'
        model.penstruct = sedumi2penbmi(interfacedata.F_struc,interfacedata.c,interfacedata.Q,interfacedata.K,interfacedata.monomtable,interfacedata.options,interfacedata.x0);

    case 'mpt'
        interfacedata.parametric_variables = find(ismember(recoverdata.used_variables,getvariables(F(find(is(F,'parametric'))))));
        interfacedata.requested_variables = [];
        model = yalmip2mpt(interfacedata);

    case 'penbmi-penopt'
        model.penstruct = sedumi2penbmi(interfacedata.F_struc,interfacedata.c,interfacedata.Q,interfacedata.K,interfacedata.monomtable,interfacedata.options,interfacedata.x0);

    case {'qpip','qpas'}
        model = yalmip2quadprog(interfacedata);
        model.options = interfacedata.options.qpip;

    otherwise
        model = [];
end
