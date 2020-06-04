function varargout = crossentropy(varargin)
% CROSSENTROPY
%
% y = CROSSENTROPY(x,y)
%
% Computes/declares cross entropy -sum(x.*log(y))
%
% See also ENTROPY, KULLBACKLEIBLER

switch class(varargin{1})
       
    case {'sdpvar','ndsdpvar'}
        
        varargin{1} = reshape(varargin{1},[],1);
        varargin{2} = reshape(varargin{2},[],1);
        
        if length(varargin{1})~=length(varargin{2})
            if length(varargin{1})==1
                varargin{1} = repmat(varargin{1},length(varargin{2}),1);
            elseif  length(varargin{2})==1
                varargin{2} = repmat(varargin{2},length(varargin{1}),1);
            else
                error('Dimension mismatch in crossentropy')
            end
        end
        
        varargout{1} = yalmip('define','crossentropy_internal',[varargin{1};varargin{2}]);
            
    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end




