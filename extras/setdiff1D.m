function [c,ia] = setdiff(a,b,flag)
%SETDIFF Set difference.
%   SETDIFF(A,B) when A and B are vectors returns the values
%   in A that are not in B.  The result will be sorted.  A and B
%   can be cell arrays of strings.
%
%   SETDIFF(A,B,'rows') when A are B are matrices with the same
%   number of columns returns the rows from A that are not in B.
%
%   [C,I] = SETDIFF(...) also returns an index vector I such that
%   C = A(I) (or C = A(I,:)).
%
%   See also UNIQUE, UNION, INTERSECT, SETXOR, ISMEMBER.

%   Copyright 1984-2003 The MathWorks, Inc.
%   $Revision: 1.5 $  $Date: 2005-02-09 22:08:06 $

%   Cell array implementation in @cell/setdiff.m

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
  
  if prod(size(a))~=numelA | prod(size(b))~=numelB
    error('MATLAB:SETDIFF:AandBvectorsOrRowsFlag', ...
          'A and B must be vectors or ''rows'' must be specified.');
  end
  
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
    
    % Convert to double arrays, which sort faster than other types.
    
    whichclass = class(a);    
    isdouble = strcmp(whichclass,'double');
    
    if ~isdouble
      a = double(a);
    end
    
    if ~strcmp(class(b),'double')
      b = double(b);
    end
    
    % Call ISMEMBER to determine list of non-matching elements of A.
    tf = ~(ismember(a,b));
    c = a(tf);
    
    % Call UNIQUE to remove duplicates from list of non-matches.
    if nargout <= 1
      c = unique(c);
    else
      [c,ndx] = unique(c);
      
      % Find indices by using TF and NDX.
      where = find(tf);
      ia = where(ndx);
    end
    
    % Re-convert to correct output data type using FEVAL.
    if ~isdouble
      c = feval(whichclass,c);
    end    
  end
  
  % If row vector, return as row vector.
  if rowvec
    c = c.';
    if nOut > 1
      ia = ia.';
    end
  end
  

% Automatically deblank strings
if ischar(a)
  c = deblank(c);
end
