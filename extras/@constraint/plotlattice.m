function varargout = plotlattice(varargin)
%PLOTLATTICE  Plots an integer lattice
%
% p = plotlattice(C,which,c,size,options)
%
% Note that only convex sets C in R^2 are supported.
%
% C    :  Constraint object
% which:  'inner' or 'outer'
% color:  color [double] ([r g b] format) or char from 'rymcgbk'
% size :  Size of marker
% options: options structure from sdpsettings

varargin{1} = lmi(varargin{1});
if nargout == 0
     plotlattice(varargin{:});
else
    varargout{1} = plotlattice(varargin{:});
end
