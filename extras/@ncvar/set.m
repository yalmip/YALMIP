function F = set(varargin)
%SET Defines a constraint (the feasible set)
%   
%    F = SET([])           Creates an empty SET-object
%    F = SET(X > Y)        Constrains X-Y to be positive semi-definite if X-Y is Hermitian,
%                          interpreted as element-wise constraint otherwise
%    F = SET(X==Y)         Element-wise equality constraint
%    F = SET(CONE(X,Y))    Second order cone constraint ||X||<Y (X column vector, Y scalar)
%
%  Constraints can also be generated using string notation (displays nicely with syntax high-lightning)
%    F = SET('X>Y')        Constrains X-Y to be positive semi-definite if X-Y is Hermitian,
%                          interpreted as element-wise constraint otherwise
%    F = SET('X==Y')       Element-wise equality constraint  
%    F = SET('||X||<Y')    Create second order cone constraint (X and Y column vectors)
%    F = SET('cone(X,Y)')  Create second order cone constraint (X and Y column vectors)
%
%  Variables can be constrained to be integer or binary
%    F = SET(INTEGER(X))
%    F = SET(BINARY(X))
%
%  Multiple constraints are obtained with overloaded plus
%    F = SET(X > 0) + SET(CONE(X(:),1)) + SET(X(1,1) == 1/2)
%
%  Double-sided constraint (and extensions) can easily be defined
%  The following two comands give equivalent problems
%    F = SET(X > 0 > Y > Z < 5 < W)
%    F = SET(X > 0) + SET(0 > Y) + SET(Y > Z) + SET(Z < 5) + set(5 < W)
%
%  A constraint can be tagged with a name or description 
%    F = SET(X > Y,'tag')  Gives the constraint a description (used in display/checkset)  
%
%  General info
%
%    The right-hand side and left-hand side can be interchanged. Supports {>,<,>=,<=,==}.
%    See FAQ for information on strict vs. non-strict inequalities.
%
%    Any valid expression built using DOUBLEs & SDPVARs can be used on both sides.
%
%    The advantage of using the string notation approach is that more information will be  
%    shown when the SET is displayed (and in checkset)
%
%    See also DUAL, SOLVESDP, INTEGER, BINARY



if isa(varargin{1},'blkvar')
    varargin{1} = sdpvar(varargin{1});     
end

switch nargin
case 0
    F = lmi;
case 1
    F = lmi(varargin{1});
case 2
    F = lmi(varargin{1},varargin{2});
case 3
    F = lmi(varargin{1},varargin{1},varargin{3});
otherwise
end
    