function varargout = softplus(varargin)
%SOFTPLUS
%
% y = SOFTPLUS(x)
%
% Computes/declares log(1 + exp(x)) ( = logsumexp([0;x]) )
%
% Implemented as evalutation based nonlinear operator. Hence, the convexity
% of this function is exploited to perform convexity analysis and rigorous
% modelling.

switch class(varargin{1})

    case 'double'
        x = varargin{1};
        varargout{1} = log(1 + exp(x));

    case 'sdpvar'

         x = varargin{1};
         if numel(x)>1           
            y = [];
            for i = 1:numel(x)
                y = [y;logsumexp([0;x(i)])];
            end            
            varargout{1} = reshape(y,size(x));          
        else
            varargout{1} = logsumexp([0;varargin{1}]);
        end
   
    otherwise
        error([upper(mfilename) ' called with weird argument']);
end