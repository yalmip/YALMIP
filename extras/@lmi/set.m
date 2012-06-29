function F = set(F,tag)
%set               Defines a constraint (the feasible set)
%   
%    F = SET([])           Creates an empty SET-object
%    F = SET(X > Y)        Constrains X-Y to be positive semi-definite if X-Y is Hermitian,
%                          interpreted as element-wise constraint otherwise
%    F = SET(X==Y)         Element-wise equality constraint
%    F = SET(CONE(X,Y))    Second order cone constraint ||X||<Y (X and Y column vectors)
%
%  Constraints can also be generated using string notation (displays nicely with syntax high-lightning)
%    F = SET('X>Y')        Constrains X-Y to be positive semi-definite if X-Y is Hermitian,
%                          interpreted as element-wise constraint otherwise
%    F = SET('X==Y')       Element-wise equality constraint  
%    F = SET('||X||<Y')    Create second order cone constraint (X and Y column vectors)
%
%  One can also use overloaded >, < and ==
%
%  Variables can be constrained to be integer or binary
%    F = SET(INTEGER(X))
%    F = SET(BINARY(X))
%
%  Multiple constraints are obtained with overloaded plus
%    F = set(X > 0) + set(CONE(X(:),1)) + SET(X(1,1) == 1/2)
%
%  Double-sided constraint (and extensions) can easily be defined
%  The following two comands give equivalent problems
%    F = set(X > 0 > Y > Z < 5 < W)
%    F = set(X > 0) + set(0 > Y) + set(Y > Z) + set(Z < 5) + set(5 < W)
%
%  A constraint can be tagged with a name or description 
%    F = SET(X > Y,'tag')  Gives the constraint a tag (used in display/checkset)  
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
%    See also   DUAL, SOLVESDP, INTEGER, BINARY
if nargin == 2
    for i = 1:length(F.clauses)
        F.clauses{i}.handle = tag;
    end
end