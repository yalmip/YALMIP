function varargout = lambertw(varargin)
%LAMBERTW (overloaded)

switch class(varargin{1})

    case 'sdpvar'

        varargout{1} = InstantiateElementWise(mfilename,varargin{:});
                
    case 'char'

        switch varargin{1}
            case 'graph'
                t = varargin{2};
                x = varargin{3};                  
                sdpvar z
                F = [expcone([z;t;x]), cone([2*t;1-z],1+z)];                
                operator = CreateBasicOperator('concave','increasing','positive','graph'); 
                operator.domain = [0 inf];
            case 'exact'
                x = varargin{3};      
                F = (x >= 0);
                operator = CreateBasicOperator('concave','increasing','positive','callback');            
                operator.derivative = @(x)(1./(exp(lambertw(x)) + x));
                operator.inverse = @(x)(x.*exp(x));
                operator.domain = [0 inf];
        end

        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = x;

    otherwise
        error('SDPVAR/LAMBERTW called with CHAR argument?');
end