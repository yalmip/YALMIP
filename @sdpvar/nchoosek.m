function varargout=nchoosek(varargin)
%NCHOOSEK (overloaded)

switch class(varargin{1})
    
    case 'double'
        error('Overloaded SDPVAR/NCHOOSEK CALLED WITH DOUBLE. Report error')
        
    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        
        if nargin ~=2
            error('Expecting two arguments')
        end
        
        n = varargin{1};
        k = varargin{2};
        if isa(k,'sdpvar')
            error('SDPVAR/NCHOOSEK currently only supported for fixed k');
        end
        
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char' % YALMIP send 'graph' when it wants the epigraph or hypograph
        
        t = varargin{2};
        n = varargin{3};
        k = varargin{4};
        
        [U,L] = derivebounds(n);
        if isinf(U) | isinf(L)
            error('Bounds required on variables in nhoosek');
        end
        L = ceil(L);
        U = floor(U);
        range = L:U;
        for i = range
            if i>=k
                nchoosekVal(i) = nchoosek(i,k);
            else
                nchoosekVal(i) = 0;
            end
        end
        nchoosekVal = nchoosekVal(L:U);
        z = binvar(length(range),1);
        F = [n == range*z, t == nchoosekVal*z,sum(z)==1];
        
        varargout{1} = F;
        varargout{2} = struct('convexity','none','monotonicity','increasing','definiteness','positive','model','integer');
        varargout{3} = n;
    otherwise
        error('Strange type on first argument in SDPVAR/SORT');
end