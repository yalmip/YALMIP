function sys = probability(c)
% PROBABILITY Create basis for chance constraint
%
% EXAMPLE:
% The following example computes the largest value t such that the
% probability that a zero mean unit variance of a Gaussian variable is
% larger than t, is larger than 0.9  
%
%  w = sdpvar(1,1);
%  t = sdpvar(1);
%  Model = [probability(a >= t) >= 0.9,uncertain(a,'normal',0,1)];
%  optimize(Model,-t)

superiorto('double');
superiorto('sdpvar');
if isa(c,'constraint') | isa(c,'lmi')
    if ~is(c,'elementwise')
        error('Probability constraints only applicable to single elementwise constraints');
    else
        sys.Weight{1} = 1;
        sys.Risk{1} = sdpvar(1);
        sys.Offset{1} = 0;
        sys.Constraint{1} = c;        
        sys = class(sys,'probability');
    end
else
    error('Probability constraints only applicable to single elementwise constraints');
end