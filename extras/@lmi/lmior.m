function varargout = or(varargin)
%OR (overloaded)

% Author Johan Löfberg 
% $Id: lmior.m,v 1.2 2007-08-02 19:33:16 joloef Exp $   

% Models OR using a nonlinear operator definition
switch class(varargin{1})
    case 'char'
        z = varargin{2};
        X = varargin{3};
        Y = varargin{4};        
      
        F = set([]);                
        x = binvar(1,1); 
        F = F + set(implies_internal(x,X));
        y = binvar(1,1); 
        F = F + set(implies_internal(y,Y));        
     
        xy = [x y];      
        varargout{1} = F + set(sum(xy) >= 1);
        varargout{2} = struct('convexity','none','monotonicity','exact','definiteness','none');
        varargout{3} = recover(depends(F));

    case {'lmi'}
        x = varargin{1};
        y = varargin{2};
        varargout{1} = set(yalmip('define','lmior',varargin{:}) == 1);

    otherwise
end

end