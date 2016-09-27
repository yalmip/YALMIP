function x = uncertain(x,varargin)
%UNCERTAIN Declares a variable as uncertain
%
%   F = UNCERTAIN(W) is used to describe the set of uncertain variables
%   in an uncertain program
%
%   INPUT
%    W : SDPVAR object or list of constraints
%
%   OUTPUT
%    F : Constraint object
%
%   EXAMPLE
%    sdpvar x w
%    F = [x + w <= 1], W = [-0.5 <= w <= 0.5];
%    optimize([F,W,uncertain(w)],-x) 
%
%    sdpvar x w
%    F = [x + w <= 1], W = [-0.5 <= w <= 0.5];
%    optimize([F,W,uncertain(W)],-x) 
%
%
%   See also SOLVESDP, ROBUSTIFY

if nargin == 1 || ((nargin == 2) && strcmpi(varargin{1},'deterministic'))
    x.typeflag = 15;
    x.extra.distribution.name = 'deterministic';
    x = lmi(x);
else
    x.typeflag = 16;
    x.extra.distribution.name = varargin{1};
    x.extra.distribution.parameters = {varargin{2:end}};
    try
        if isa(x.extra.distribution.name,'function_handle')
            temp = {x.extra.distribution.name,x.extra.distribution.parameters{:},x.dim};            
            temp = feval(temp{:});
        else
            temp = {@random,x.extra.distribution.name,x.extra.distribution.parameters{:},x.dim};
            temp = feval(temp{:});
        end
    catch
        lasterr
        error('Cannot evaluate the uncertainty description supplied.')
    end
    x = lmi(x);                 
end