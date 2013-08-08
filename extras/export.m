function [model,recoverdata,diagnostic,interfacedata] = export(varargin)
%EXPORT  Exports YALMIP problem to solver specific format
%
%   [MODEL,RECOVERYMODEL,DIAGNOSTIC,INTERNAL] = EXPORT(F,h,options) is the common command to
%   export optimization problems of the following kind
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

% Arrrgh, new format with logdet much better, but we have to
% take care of old code, requires some testing...
varargin = combatible({varargin{:}});
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
        F = set(F);
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
           F = set([]);
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
        
        % Hack to handle CPLEX slow treatment of paramters. Remove all
        % default settings, so optimizer runs fast
        try
            o1 = cplexoptimset('cplex');
            o2 = interfacedata.options.cplex;
            n = fieldnames(o1);
            for i = 1:length(n)
                if isequal(o1.(n{i}),o2.(n{i}))
                    interfacedata.options.cplex = rmfield(interfacedata.options.cplex,n{i});                    
                end
            end
            model = yalmip2cplex(interfacedata);
        catch
        end

    case 'ecos'
        model = yalmip2ecos(interfacedata);      
        
        
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
                    
    case 'quadprog'
        model = yalmip2quadprog(interfacedata);
        
    case {'sedumi-1.05','sedumi-1.1','sedumi-1.3'}        
        model = yalmip2sedumi(interfacedata);
        
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
                
    case 'mpt'
        interfacedata.parametric_variables = find(ismember(recoverdata.used_variables,getvariables(F(find(is(F,'parametric'))))));
        interfacedata.requested_variables = [];
        model = yalmip2mpt(interfacedata);
                
    case 'penbmi-penopt'        
        model.penstruct = sedumi2penbmi(interfacedata.F_struc,interfacedata.c,interfacedata.Q,interfacedata.K,interfacedata.monomtable,interfacedata.options,interfacedata.x0);
               
    otherwise
        model = [];
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
elseif m>=3 & isstruct(varargin{3})
    classification = 2;
elseif m>=4 & isstruct(varargin{4})
    classification = 1;
elseif m>=2 & isa(varargin{2},'lmi')
    classification = 1;
elseif m>=3 & isa(varargin{3},'sdpvar')
    classification = 1;
elseif m>=2 & isa(varargin{2},'sdpvar') & min(size(varargin{2}))==1
    classification = 2;
elseif m>=2 & isa(varargin{2},'sdpvar') & prod(size(varargin{2}))>=1
    classification = 1;
elseif m>=2 & isa(varargin{2},'logdet')
    classification = 2;
elseif m==2 & isempty(varargin{2})
    classification = 2;
elseif m>=3 & isempty(varargin{2}) & isempty(varargin{3})
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



function model = yalmip2gurobimex(interfacedata)
% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
Q       = interfacedata.Q;
c       = interfacedata.c;
K       = interfacedata.K;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
semicont_variables = interfacedata.semicont_variables;
ub      = interfacedata.ub;
lb      = interfacedata.lb;
x0      = interfacedata.x0;
interfacedata.gettime = 0;
n = length(c);

if ~isempty(ub)
    LB = lb;
    UB = ub;
    LB(binary_variables)  = round(LB(binary_variables));
    UB(binary_variables)  = round(UB(binary_variables));
    LB(integer_variables) = round(LB(integer_variables));
    UB(integer_variables) = round(UB(integer_variables));
else
    LB = [];
    UB = [];
end

if options.showprogress;showprogress('Calling GUROBI',options.showprogress);end

if ~isempty(semicont_variables)
    % Bounds must be placed in LB/UB
    [LB,UB,cand_rows_eq,cand_rows_lp] = findulb(F_struc,K,LB,UB);
    F_struc(K.f+cand_rows_lp,:)=[];
    F_struc(cand_rows_eq,:)=[];
    K.l = K.l-length(cand_rows_lp);
    K.f = K.f-length(cand_rows_eq);
    
    redundant = find(LB<=0 & UB>=0);
    semicont_variables = setdiff(semicont_variables,redundant);
        
end

SENSE = 1;     % Minimize
C = full(c);   % Must be full
if size(F_struc,1)>0
    B = full(F_struc(:,1));         % Must be full
    A =-F_struc(:,2:end);
else
    B = [];
    A = [];
end

% Optimized code, make a lot of difference when you make this call 10000
% times in a branch and bound setting...
CTYPE = [char(ones(K.f,1)*61); char(ones(K.l,1)*60)];
VARTYPE = char(ones(length(c),1)*67);
VARTYPE(setdiff(integer_variables,semicont_variables)) = 'I';
VARTYPE(binary_variables)  = 'B';  % Should not happen except from bmibnb
VARTYPE(setdiff(semicont_variables,integer_variables)) = 'S';  % Should not happen except from bmibnb
VARTYPE(intersect(semicont_variables,integer_variables)) = 'N';

% Gurobi assumes semi-continuous variables only can take negative values so
% we negate semi-continuous violating this
NegativeSemiVar = [];
if ~isempty(semicont_variables)
    NegativeSemiVar = find(UB(semicont_variables) < 0);
    if ~isempty(NegativeSemiVar)
        temp = UB(semicont_variables(NegativeSemiVar));
        UB(semicont_variables(NegativeSemiVar)) = -LB(semicont_variables(NegativeSemiVar));
        LB(semicont_variables(NegativeSemiVar)) = -temp;
        A(:,semicont_variables(NegativeSemiVar)) = -A(:,semicont_variables(NegativeSemiVar));
        C(semicont_variables(NegativeSemiVar)) = -C(semicont_variables(NegativeSemiVar));
        if ~isempty(x0)
            x0(NegativeSemiVar) = -NegativeSemiVar;
        end
    end
end

if nnz(Q)>0
    [ii,jj,kk] = find(Q);
    ii = ii-1;
    jj = jj-1;
    comp = computer;
    % According to Wotao's testing
    if isempty(strfind(comp,'64')) | isequal(comp,'GLNXA64') | isequal(comp,'PCWIN64')
        options.gurobi.QP.qrow = int32(ii)';
        options.gurobi.QP.qcol = int32(jj)';
    elseif strfind(comp,'MACI64')
        options.gurobi.QP.qrow = int32(ii)';
        options.gurobi.QP.qcol = int32(jj)';    
    else
        options.gurobi.QP.qrow = int64(ii)';
        options.gurobi.QP.qcol = int64(jj)';
    end
    options.gurobi.QP.qval = kk';
end
        
if ~options.verbose
    options.gurobi.DisplayInterval = 0;
    options.gurobi.Display = 0;
end

if ~isempty(x0)
    options.gurobi.Start = x0(:)';
end

if options.savedebug
    save gurobidebug
end

if ~isempty(K.sos.type)
    options.gurobi.SOS.weights = spalloc(length(c),length(K.sos.type),0);
    for i = 1:length(K.sos.type)
        options.gurobi.SOS.types(i)= int32(str2num(K.sos.type(i)));
        options.gurobi.SOS.weights(K.sos.variables{i},i) = 1:length(K.sos.variables{i});
    end
end
model.C = C;


% Call mex-interface
solvertime = clock; 
if isempty(binary_variables) & isempty(integer_variables) & isempty(semicont_variables)
    [x,val,flag,output,lambda] = gurobi_mex(C,SENSE,sparse(A),B,CTYPE,LB,UB,VARTYPE,options.gurobi);
else
    [x,val,flag,output] = gurobi_mex(C,SENSE,sparse(A),B,CTYPE,LB,UB,VARTYPE,options.gurobi);
    lambda = [];
end
