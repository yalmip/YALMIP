function Y = max(X,aux,workdim)
% abs (overloaded)

switch nargin
    case 1
        workdim = 1;
    case 2
        Y = max(reshape(X,[],1),reshape(aux,[],1));
        Y = reshape(Y,size(X));
        return
    case 3
        if ~isempty(aux)
            error('MAX with two matrices to compare and a working dimension is not supported in MATLAB.');
        end
    otherwise
        error('Too many input arguments')
end

xdim = size(X);

if workdim > length(xdim)
    error('Index exceeds matrix dimensions.');
end

newdim = xdim;
newdim = circshift(xdim',-(workdim-1))';
newdim(1)=1;
Y = reshape(max(reshape(shiftdim(X,workdim-1),xdim(workdim),[])),newdim);
Y = shiftdim(Y,length(xdim)-workdim+1);
