function X = setcutflag(X,flag)
%setcutlag Intenal : defines a SET object as a CUT
%
% Author Johan Löfberg
% $Id: setcutflag.m,v 1.3 2005-02-10 12:26:38 johanl Exp $
if nargin == 1
    flag = 1;
end
X.clauses{1}.cut = flag;