function monom = genmonoms(newton_m,x);
%GENMONOMS Internal function used in SOS programs

% Author Johan Löfberg
% $Id: genmonoms.m,v 1.1 2006-03-30 13:56:54 joloef Exp $


precalc = [];
monom = [];
for i = 1:size(newton_m,2)
    powers = unique(newton_m(:,i));
    for j = 1:length(powers)         
        precalc{1+powers(j),i} = x(i)^powers(j);
    end   
end

monom = [];
for i = 1:size(newton_m,1)
    temp = 1;
    for j = 1:size(newton_m,2)
        temp = temp*precalc{newton_m(i,j)+1,j};%x(j)^newton_m(i,j);
    end
    monom = [monom;temp];
end