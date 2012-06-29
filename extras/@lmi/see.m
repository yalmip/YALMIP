function see(F)
%see               Displays internal structure of matrix variable in a constraint

% Author Johan Löfberg
% $Id: see.m,v 1.3 2005-02-08 16:11:17 johanl Exp $

if length(F.clauses)>1
    disp('Can only apply SEE on set objects with one constraint');
else
    see(F.clauses{1}.data);
end