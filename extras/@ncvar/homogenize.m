function p_homo = homogenize(p,y);
%HOMOGENIZE Homogenize polynomial
%
% f = homogenize(p,x)

% Author Johan Löfberg 
% $Id: homogenize.m,v 1.1 2006-08-10 18:00:20 joloef Exp $  

deg   = degree(p);
deg_y = degree(y);
if rem(deg,deg_y)~=0
    error('The degree of the homogenizer is not an even fraction of deg(p).');
end

if 0
    error('The homogenizer must be homogenious.');
end

p_variables = getvariables(p);
p_homo = getbasematrix(p,0)*y^(deg/deg_y);
for i = 1:length(p_variables);
    monom = recover(p_variables(i));
    if degree(monom)<deg
        power = (deg-(degree(monom)))/deg_y;
        p_homo = p_homo + getbasematrix(p,p_variables(i))*monom*y^power;
    else
       p_homo = p_homo + getbasematrix(p,p_variables(i))*monom;
   end;
end
% Reset info about conic terms
p_homo.conicinfo = [0 0];  