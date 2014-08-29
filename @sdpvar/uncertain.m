function x = uncertain(x,varargin)
%UNCERTAIN Declares a variable as uncertain
%
%   F = UNCERTAIN(w) is used to describe the set of uncertain variables
%   in an uncertain program
%
%   INPUT
%    w : SDPVAR object
%
%   OUTPUT
%    F : Constraint object
%
%   EXAMPLE
%    sdpvar x w
%    F = [x + w <= 1, -0.5 <= w <= 0.5, uncertain(w)];
%    solvesdp(F,-x) 
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
    x = lmi(x);
end