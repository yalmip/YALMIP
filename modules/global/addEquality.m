function p = addEquality(p,row)
% Equalities on top so easy
p.F_struc = [row;p.F_struc];
p.K.f = p.K.f + size(row,1);