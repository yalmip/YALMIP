function varargout = exp(varargin)

switch class(varargin{1})

    case 'sdpvar'
        x = varargin{1};
        d = size(x);
        x = x(:);
        y = [];
        for i = 1:prod(d)
            xi = extsubsref(x,i);
            if isnumeric(xi)
                y = [y;exp(xi)];
            else
                if isreal(xi)
                    if i>1
                        y = [y;InstantiateElementWise(mfilename,xi)];
                    else
                        y = InstantiateElementWise(mfilename,xi);
                    end
                else
                    re = real(xi);
                    im = imag(xi);
                    y = [y;exp(re)*(cos(im) + sqrt(-1)*sin(im))];
                end
            end
        end
        varargout{1} = reshape(y,d);
                    
    case 'char'
        
        operator = CreateBasicOperator('convex','increasing','positive','callback');
        operator.derivative = @(x)exp(x);
        operator.inverse = @(x)log(x);   
        operator.range = [0 inf];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end