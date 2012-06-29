function y = subsref(X,Y)
%SUBSREF (overloaded)

% Author Johan Löfberg 
% $Id: subsref.m,v 1.2 2005-10-13 14:00:43 joloef Exp $   

if length(Y.subs)>=1
    if isequal(Y.subs{1},':')
        Y.subs{1} = 1:X.n;
    end
end
if length(Y.subs)>=2
    if isequal(Y.subs{2},':')
        Y.subs{2} = 1:X.m;
    end
end

if length(Y.subs{1}) < length(Y.subs{2})
    
    [usedCols,vals] = ismember(X.jX,Y.subs{2});
    X.jX = vals(usedCols);
    X.iX = X.iX(usedCols);
    X.sX = X.sX(usedCols);
       
    [usedRows,vals] = ismember(X.iX,Y.subs{1});
    X.iX = vals(usedRows);
    X.jX = X.jX(usedRows);
    X.sX = X.sX(usedRows);    
    
else
   
    [usedRows,vals] = ismember(X.iX,Y.subs{1});
    X.iX = vals(usedRows);
    X.jX = X.jX(usedRows);
    X.sX = X.sX(usedRows);
    
    [usedCols,vals] = ismember(X.jX,Y.subs{2});
    X.jX = vals(usedCols);
    X.iX = X.iX(usedCols);
    X.sX = X.sX(usedCols);
    
end
    
y = sparse(X.iX,X.jX,X.sX,length(unique(Y.subs{1})),length(unique(Y.subs{2})));


function [col,row] = getXYfromi(n,s)

row = ceil(min(roots([-1/2 n+1/2 -s])));

if row == 1
    col = s;
else
   row = ceil(-1e-6+min(roots([-1/2 n+1/2 -s])));
    i = row-1;
    i = (i)*n-((i-1)*(i)/2) + 1;
    col = (s-i)+row;
end
    

