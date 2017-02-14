function varargout = plot(varargin)
%PLOT  Plots the feasible region of a set of constraints
%
% p = plot(C,x,c,n,options)
%
% Note that only convex sets are allowed, or union of convex sets
% represented using binary variables (either defined explictly or
% introduced by YALMIP when modelling, e.g., mixed integer linear
% programming representable operators)
%
% C:  Constraint object
% x:  Plot variables [At most three variables]
% c:  color [double] ([r g b] format) or char from 'rymcgbk'
% n:  #vertices [double ] (default 100 in 2D and 300 otherwise)
% options: options structure from sdpsettings

% Get the onstraints
if nargin<1
    return
end

F = varargin{1};

if length(F)==0
    return;
end

if nargin < 5
    opts = sdpsettings('verbose',0);
else
    opts = varargin{5};
    if isempty(opts)
        opts = sdpsettings('verbose',0);
    end
end
opts.verbose = max(opts.verbose-1,0);

if any(is(F,'uncertain'))
    F = robustify(F,[],opts);
end

if nargin < 3
    color=['rymcgbk']';
else
    color = varargin{3};
    if isa(color,'sdpvar')
        error('The variables should be specified in the second argument.');
    end
    color = color(:)';
    if isempty(color)
        color=['rymcgbk']';
    end
end

% Plot onto this projection (at most in 3D)
if nargin < 2
    x = [];
else
    x = varargin{2};
    if ~isempty(x)
        if ~isa(x,'sdpvar')
            error('Second argument should be empty or an SDPVAR');
        end
        x = x(:);
        x = x(1:min(3,length(x)));
    end    
end

if isempty(F)
    return
end

if any(is(F,'sos'))
    % Use image representation, feels safer (I'm not sure about the logic
    % in the code at the moment. Can the dualization mess up something
    % otherwise...)
    if ~(opts.sos.model == 1)
        opts.sos.model = 2;
    end
    % Assume the variables we are plotting are the parametric
    F = compilesos(F,[],opts,x);
end

% Create a model in YALMIPs low level format
% All we change later is the cost vector
%sol = solvesdp(F,sum(x),opts);
[model,recoverdata,diagnostic,internalmodel] = export(F,[],opts,[],[],0);
if isempty(internalmodel) | (~isempty(diagnostic) && diagnostic.problem)
    error('Could not create model. Can you actually solve problems with this model?')
end
internalmodel.options.saveduals = 0;
internalmodel.getsolvertime = 0;
internalmodel.options.dimacs = 0;

% Try to find a suitable set to plot
if isempty(x)
    if isempty(internalmodel.extended_variables) & isempty(internalmodel.aux_variables)
        x = depends(F);
        x = x(1:min(3,length(x)));
        localindex = 1;
        localindex = find(ismember(recoverdata.used_variables,x));
    else
        % not extended variables
        candidates = setdiff(1:length(internalmodel.c),[ internalmodel.aux_variables(:)' internalmodel.extended_variables(:)']);
        % Not nonlinear variables
        candidates = candidates(find(internalmodel.variabletype(candidates)==0));
        % Settle with this guess
        localindex = candidates(1:min(3,length(candidates)));
        x = localindex;
    end
else
    localindex = [];
    for i = 1:length(x)
        localindex = [localindex find(ismember(recoverdata.used_variables,getvariables(x(i))))];
    end
end

if nargin < 4
    if length(x)==3
        n = 300;
    else
        n = 100;
    end
else
    n = varargin{4};
    if isempty(n)
        if length(x)==3
            n = 300;
        else
            n = 100;
        end
    end
    if ~isa(n,'double')
        error('4th argument should be an integer>0');
    end
end

if ~isempty(internalmodel.integer_variables)
    error('PLOT can currently not display sets involving integer variables');
end

% Use common function for lmi/plot and optimizer/plot
if nargout == 0
    plotInternalModel(internalmodel,x,n,localindex,color,opts);
else
    varargout{1} = plotInternalModel(internalmodel,x,n,localindex,color,opts);
end