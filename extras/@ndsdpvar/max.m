function Y = max2(X,aux,workdim)
% abs (overloaded)
xdim = size(X);

switch nargin
    case 1
        sdim = 1;
        workdim=find(xdim>1,1);

    case 2
        Y = max(reshape(X,[],1),reshape(aux,[],1));
        Y = reshape(Y,xdim);
        return
    case 3
        
        if ~isempty(aux)
            error('MAX with two matrices to compare and a working dimension is not supported in MATLAB.');
        elseif xdim(workdim)==1
            Y=X;
            return
        else
            sdim=workdim;
        end
        
    otherwise
        error('Too many input arguments')
end


if workdim > length(xdim)
    error('Index exceeds matrix dimensions.');
end

newdim = xdim;
newdim(workdim)=1;
newdim = circshift(newdim',-(sdim-1))';

Y = reshape(max(reshape(shiftdim(X,workdim-1),xdim(workdim),[])),newdim);
Y = shiftdim(Y,length(xdim)-sdim+1);
