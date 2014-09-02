function varargout = plot(varargin)
%PLOT  Plots the feasible region of a set of constraints
%
% p = plot(C,x,c,n,options)
%
% Note that only convex sets are allowed, or union of convex sets
% represented using binary variables (either defined explictly or
% introduced by YALMIP when modelling, e.g., mixed integer linear
% programming represetable operators)
%
% C:  Constraint object
% x:  Plot variables [At most three variables]
% c:  color [double] ([r g b] format) or char from 'rymcgbk'
% n:  #vertices [double ]
% options: options structure from sdpsettings

varargin{1} = lmi(varargin{1});
if nargout == 0
     plot(varargin{:});
else
    varargout{1} = plot(varargin{:});
end
