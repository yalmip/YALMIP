function varargout=quadratic_over_affine(varargin)
%quadratic_over_affine (overloaded)
%
% p = quadratic_over_affine(x,t)
%
% Returns p = (x.^2)./t

switch class(varargin{1})
    case 'double'
        varargout{1} = varargin{1}/varargin{2};
        
    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        if is(varargin{1},'quadratic') & ~is(varargin{1},'linear') & is(varargin{2},'linear')
            p = varargin{1};
            t = varargin{2};
            if length(p)>1
                if length(t)==1
                    t = repmat(t,size(p,1),size(p,2));
                end
                temp = [];
                for i = 1:length(p)
                    temp = [temp;quadratic_over_affine(p(i),t(i))];
                end
                temp = reshape(temp,size(p,1),size(p,2));
                varargout{1} = temp;
            elseif length(t)>1
                temp = [];
                for i = 1:length(t)
                    temp = [temp;quadratic_over_affine(p,t(i))];
                end
                temp = reshape(temp,size(t,1),size(t,2));
                varargout{1} = temp;
            else
                [Q,c,f,x,info] = quaddecomp(p);
                q = chol(Q)*x;
                varargin{1} = q;
                varargout{1} = yalmip('define','quadratic_over_affine_expanded',varargin{:});
            end
        else
            error('First argument should be quadratic and second affine')
        end
        
        
    otherwise
        error([upper(mfilename) ' called with weird argument']);
end
