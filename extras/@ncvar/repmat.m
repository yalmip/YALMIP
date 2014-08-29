function Y=repmat(varargin)
%REPMAT (overloaded)

try
  X = varargin{1};
  Y = X;
  Y.basis = [];
  n = Y.dim(1);
  m = Y.dim(2);
  for i = 1:length(Y.lmi_variables)+1
    temp = repmatfixed(reshape(X.basis(:,i),n,m),varargin{2:end});
    Y.basis(:,i) = temp(:);
  end
  Y.dim(1) = size(temp,1);
  Y.dim(2) = size(temp,2);
  % Reset info about conic terms
  Y.conicinfo = [0 0];
catch
  error(lasterr)
end



function B = repmatfixed(A,M,N)

if nargin < 2
    error('MATLAB:repmat:NotEnoughInputs', 'Requires at least 2 inputs.')
end

if nargin == 2
    if isscalar(M)
        siz = [M M];
    else
        siz = M;
    end
else
    siz = [M N];
end

if isscalar(A)
    nelems = prod(siz);
    if nelems>0
        % Since B doesn't exist, the first statement creates a B with
        % the right size and type.  Then use scalar expansion to
        % fill the array. Finally reshape to the specified size.
        B = spalloc(nelems,1,nnz(A));
        B(nelems) = A;
        if ~isequal(B(1), B(nelems)) | ~(isnumeric(A) | islogical(A))
            % if B(1) is the same as B(nelems), then the default value filled in for
            % B(1:end-1) is already A, so we don't need to waste time redoing
            % this operation. (This optimizes the case that A is a scalar zero of
            % some class.)
            B(:) = A;
        end
        B = reshape(B,siz);
    else
        B = A(ones(siz));
    end
elseif ndims(A) == 2 & numel(siz) == 2
    [m,n] = size(A);
    
    if (m == 1 & siz(2) == 1)
        B = A(ones(siz(1), 1), :);
    elseif (n == 1 & siz(1) == 1)
        B = A(:, ones(siz(2), 1));
    else
        mind = (1:m)';
        nind = (1:n)';
        mind = mind(:,ones(1,siz(1)));
        nind = nind(:,ones(1,siz(2)));
        B = A(mind,nind);
    end
else
    Asiz = size(A);
    Asiz = [Asiz ones(1,length(siz)-length(Asiz))];
    siz = [siz ones(1,length(Asiz)-length(siz))];
    for i=length(Asiz):-1:1
        ind = (1:Asiz(i))';
        subs{i} = ind(:,ones(1,siz(i)));
    end
    B = A(subs{:});
end

function a = isscalar(b)
[n,m] = size(b);
a = (n*m == 1);