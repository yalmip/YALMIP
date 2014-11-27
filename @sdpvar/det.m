function varargout = det(varargin)
%DET (overloaded)
%
% t = DET(X)

switch class(varargin{1})
       
    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        
        X = varargin{1};
        [n,m] = size(X);
        if n~=m
            error('Matrix must be square.')
        end
        if nargin == 2
            if strcmp(varargin{2},'polynomial')
                varargout{1} = polynomialform(X);
                return
            else
                error('If you use two arguments in @sdpvar/det, the second should be ''polynomial''');
            end
        end
        if n==1
            varargout{1} = X;
            return
        else
            y = yalmip('define','det_internal',reshape(X,[],1));
        end
        varargout{1} = y;

    otherwise
end


function d = polynomialform(X)

n = X.dim(1);
m = X.dim(2);

if n~=m
    error('Matrix must be square.');
else
    switch n
        case 1
            d = X;
        case 2
            % Freakin overloading on multiplication doesn't work. Probalby
            % stupid code...
            Y1.type = '()';
            Y2.type = '()';
            Y3.type = '()';
            Y4.type = '()';
            Y1.subs = {1,1};
            Y2.subs = {2,2};
            Y3.subs = {1,2};
            Y4.subs = {2,1};
            d = subsref(X,Y1)*subsref(X,Y2)-subsref(X,Y3)*subsref(X,Y4);
        otherwise
            d = 0;
            Y.type = '()';
            for i = 1:n
                Y.subs = {i,1};
                xi = subsref(X,Y);
                if ~isequal(xi,0)
                    Y.subs = {[1:1:i-1 i+1:1:n],2:n};
                    subX = subsref(X,Y);
                    if isa(subX,'sdpvar')
                        d = d + (-1)^(i+1)*xi*polynomialform(subX);
                    else
                        d = d + (-1)^(i+1)*xi*det(subX);
                    end
                end
            end
    end
end
% Reset info about conic terms
if isa(d,'sdpvar')
    d.conicinfo = [0 0];
end