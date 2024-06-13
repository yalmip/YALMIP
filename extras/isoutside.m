function varargout = isoutside(varargin)
%ISOUTSIDE Define avoidance constraint on polytope
%
% F = ISOUTSIDE(P)
%
% Input
%  P : Linear inequalities defining polytope
%
% Since the modelling is based on big-M formulations, all
% involved variables should be explicitly bounded.

% Here is the real overloaded ismember
switch class(varargin{1})
    case {'lmi','constraint'}
        
        P = varargin{1};
        if ~(all(is(P,'linear')) & all(is(P,'elementwise')))
            error('ISOUTSIDE only applicable to linear inequalities')
        end
        varargout{1} = (yalmip('define',mfilename,varargin{:}) == 1);
        varargout{1} = setupMeta(lmi([]), mfilename,varargin{:});
        
    case 'char'
        varargout{1} = isoutside_internal(varargin{3});
end


