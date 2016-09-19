function varargout = sosmodel(varargin)
%SOSMODEL Derive sum-of-squares model without solving it
%
%    [F,obj,m] = sosmodel(F,h,options,params,monomials) compiles the SOS
%    problem (i.e., derives the SDP model) without actually solving it
%
%    Inputs
%     F         : The model involving SOS constraints
%     h         : Objective function (function of params) [optional]
%     options   : SDPSETTINGS structure [optional]
%     params    : Parametric variables in model [optional]
%     monomials : Prespecified monomials to be used [optional]
%
%    Outputs
%     F         : Constraints defining the problem
%     h         : Objective function
%     m         : Monomials used in the decomposition
%
%
% NOTE: If you do not specify any options structure or specify sos.model
% as default 0 an image model representation will be used (sos.model = 2).
% With this, the compiled model is expressed in terms of the original
% variables (with kernel model, a dual problem is derived, and thus no
% original variables are used)
% 
% See also SOLVESOS, OPTIMIZE, SOS, OPTIMIZER

% When sosmodel is used, default is to return a model expressed in the
% original variables, i.e., image space model
if nargin > 3
    ops = varargin{3};
    if isempty(ops)
        ops = sdpsettings('sos.model',2);
        varargin{3} = ops;
    else
        try
            if isequal(ops.sos.model,0)
                ops.sos.model = 2;
            end
            varargin{3} = ops;
        catch
        end
    end
else
    ops = sdpsettings('sos.model',2);
    varargin{3} = ops;
end
[varargout{1:nargout}] = compilesos(varargin{:});