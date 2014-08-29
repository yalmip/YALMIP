function [c] = setdiff1D(a,b)

nIn = 2;
isrows = 0;
rowsA = size(a,1);
colsA = size(a,2);
rowsB = size(b,1);
colsB = size(b,2);
rowvec = ~((rowsA > 1 & colsB <= 1) | (rowsB > 1 & colsA <= 1) | isrows);
nOut = 1;
numelA = length(a);
numelB = length(b);
% Handle empty arrays.
if (numelA == 0)
    % Predefine outputs to be of the correct type.
    c = a([]);
    ia = [];
    % Ambiguous if no way to determine whether to return a row or column.
    ambiguous = (rowsA==0 & colsA==0) & ...
        ((rowsB==0 & colsB==0) | numelB == 1);
    if ~ambiguous
        c = reshape(c,0,1);
        ia = reshape(ia,0,1);
    end
elseif (numelB == 0)
    % If B is empty, invoke UNIQUE to remove duplicates from A.
    if nOut <= 1
        c = unique(a);
    else
        [c,ia] = unique(a);
    end
    return
    
    % Handle scalar: one element.  Scalar A done only.
    % Scalar B handled within ISMEMBER and general implementation.
    
elseif (numelA == 1)
    if ~ismember(a,b)
        c = a;
        ia = 1;
    else
        c = [];
        ia = [];
    end
    return
    
    % General handling.
    
else
    
    % Convert to columns.
    a = a(:);
    b = b(:);
    %a = double(a);
    %b = double(b);
    %  end
    
    % Call ISMEMBER to determine list of non-matching elements of A.
    tf = ~(ismember(a,b));
    %tf = ~(ismembc(sort(a),sort(b)));
    c = a(tf);
    c = uniquestripped(c);
    
end

% If row vector, return as row vector.
if rowvec
    c = c.';
    if nOut > 1
        ia = ia.';
    end
end

