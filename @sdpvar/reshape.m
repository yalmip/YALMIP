function Y=reshape(varargin)
%RESHAPE (overloaded)

try
    Y = varargin{1};
    n = length(Y.lmi_variables);
    
    if (length(varargin{2}) <= 2) && nargin <=3
        % LAZY....
        if (Y.dim(1)==1) && (Y.dim(2)==1)
            % BUG IN MATLAB R13 ON SPARSE SCALAR IN RESHAPE!!!!!
            temp = reshape(reshape(full(Y.basis(:,1)),Y.dim(1),Y.dim(2)),varargin{2:end});
        else
            temp = reshape(reshape(Y.basis(:,1),Y.dim(1),Y.dim(2)),varargin{2:end});
        end

        % RESHAPE DOES NOT CHANGE INTERNAL DATA
        Y.dim(1) = size(temp,1);
        Y.dim(2) = size(temp,2);
        % Reset info about conic terms
        Y.conicinfo = [0 0];        
    else
        Y = reshape(ndsdpvar(Y),varargin{2:end});
    end
catch
    error(lasterr)
end
