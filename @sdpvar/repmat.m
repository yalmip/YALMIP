function Y=repmat(varargin)
%REPMAT (overloaded)

try
  X = varargin{1};
  Y = X;
  Y.basis = [];
  dim = Y.dim;
  for i = 1:length(Y.lmi_variables)+1
    temp = repmatfixed(reshape(X.basis(:,i),dim),varargin{2:end});
    Y.basis(:,i) = temp(:);
  end
  Y.dim = size(temp); 
  Y = flush(Y);
  % Reset info about conic terms
  Y.conicinfo = [0 0];
  if length(Y.dim)>2     
      Y = ndsdpvar(Y);
  end
catch
  error
end

function B = repmatfixed(A,siz,sz2)

if nargin == 3
    siz = [siz sz2];
end

% nd-indexing does not work on sparse
if length(siz)>2
    A = full(A);
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