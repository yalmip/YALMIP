function [y,loc] = matrix_sdpvar_sort(varargin)

% Called from sort to do matrix sorts using repeated vector sorts
X = varargin{1};
if nargin > 1
    dim = varargin{2};
else
    dim = 1;
end

if dim == 2
    [y,loc] = matrix_sdpvar_sort(X',1);
    y = y';
    loc = loc';
else
    y = [];
    loc = [];
    for i = 1:size(X,2)
        [yi,loci] = sort(X(:,i));
        y = [y yi];
        loc = [loc loci];
    end
end
