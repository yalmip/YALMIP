function x = uncertain(x,varargin)
%UNCERTAIN Declares a variable as uncertain
%
%   F = UNCERTAIN(W) is used to describe the set of uncertain variables
%   in an uncertain program
%
%   INPUT
%    W : SDPVAR object or list of constraints
%    S : Optional distribution information for random uncertainty
%
%   OUTPUT
%    F : Constraint object
%
%   Uncertain is used to declare uncertain variables in robust
%   deterministic worst-case optimization. It can also be used to specify
%   uncertain random variables, and their associated distribution, to be
%   used in OPTIMIZER objects with the SAMPLE command.
%
%   EXAMPLE
%    
%    Robust worst-case optimization
%
%    sdpvar x w
%    F = [x + w <= 1], W = [-0.5 <= w <= 0.5];
%    optimize([F,W,uncertain(w)],-x) 
%
%    sdpvar x w
%    F = [x + w <= 1], W = [-0.5 <= w <= 0.5];
%    optimize([F,uncertain(W)],-x) 
%
%   To specify random uncertainties, you specify the distribution, and all
%   distribution parameters following the syntax n the RANDOM command in
%   the Statistics Toolbox
%
%    sdpvar x w
%    F = [x + w <= 1, uncertain(w, 'uniform',0,1)];
%    P = optimizer([F,W,uncertain(w)],-x,[],w,x)
%    S = sample(P,10); % Sample ten instances and concatenate models
%    S([])             % Solve and return optimal x
%
%    Alternatively, you can specify a function handle which generates
%    samples. YALMIP will always send a trailing argument with dimensions
%
%    F = [x + w <= 1, uncertain(w,@mysampler,myarguments1,...)];
%
%    The standard random case above would thus be recovered with
%
%    F = [x + w <= 1, uncertain(w,@random,0,1)];
%
%
%   See also OPTIMIZE, ROBUSTMODEL, OPTIMIZER, SAMPLE

if nargin == 1 || ((nargin == 2) && strcmpi(varargin{1},'deterministic'))
    x.typeflag = 15;
    x.extra.distribution.name = 'deterministic';
    x = lmi(x);
else
    x.typeflag = 16;    
    if isa(varargin{1},'function_handle')
         temp = {varargin{:},x.dim};
    else
        temp = {@random,varargin{:},x.dim};
    end
    x.extra.distribution.name = temp{1};
    x.extra.distribution.parameters = {temp{2:end-1}};
    try
        temp = feval(temp{:});        
    catch
        disp(lasterr);
        error('Trial evaluation of attached sample generator failed.')
    end
    x = lmi(x);                 
end