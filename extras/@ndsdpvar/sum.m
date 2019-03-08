function X = sum(varargin)
% SUM (overloaded)

Y = varargin{1};
X = Y;
%X.basis = [];
if nargin == 2 && isequal(varargin{2},length(X.dim))
    % smash slices
    X.basis = kron(ones(1,X.dim(end)),speye(prod(X.dim(1:end-1))))*Y.basis;
    X.dim = X.dim(1:end-1);
else
    if nargin == 1
        index = find(X.dim~=1);
        if isempty(index)
            index = length(X.dim);
        else
            index = index(1);
        end
    else
        index = varargin{2};
        if length(index) > 1
            error('Dimension argument must be a positive integer scalar within indexing range.');
        end
    end
    if index > length(X.dim)
        return
    end
    % Permute to the case which we can do fast
    i = 1:length(X.dim);
    p = circshift(i',length(X.dim)-(index))';
    X = permute(X,p);
    % Permute might have squeezed to 2 array
    if length(size(X))~=length(size(Y))
        % Expand back last dimension
        if length(size(X))==2
            X = ndsdpvar(X);
        end
        X.dim = [X.dim ones(1,length(size(Y))-length(size(X)))];
    end
    X = sum(X,length(X.dim));
    if isa(X,'sdpvar')
        X = ndsdpvar(X);
        X.dim = [X.dim ones(1,length(size(varargin{1}))-length(size(X)))];
    else
        X.dim = [X.dim ones(1,length(size(varargin{1}))-length(size(X)))];
    end
    p = circshift((1:length(X.dim))',-(length(Y.dim)-(index)))';
    X = permute(X,p);
   % X.dim = Y.dim;
   % X.dim(index) = 1;
end
if isa(X,'ndsdpvar')
    X.conicinfo = [0 0];
end
X = clean(X);