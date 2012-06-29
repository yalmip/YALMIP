function [w,F,DefinedMonoms] = recdef(pow,F,DefinedMonoms,setinitials);

% Author Johan Löfberg
% $Id: recdef.m,v 1.4 2006-03-08 16:12:51 joloef Exp $

% Recursively define monomial x(pow(1))*x(pow(2))*x(pow(3))..
% using a set of linear variables and bilinear constraints
switch length(pow)
    case 1
        % Just a linear variable, no constraints needed
        w = recover(pow);
    case 2
        % Quadratic
        w = sdpvar(1,1);
        xx = prod(recover(pow));
        if setinitials
            assign(w,double(xx));%prod(recover(pow))));
        end
        F = F + set(xx == w);        
    otherwise
        % FIX: special case x^3, implemented just to
        % speed up a particular simulation I had to run.
        % This whole file should be revised, booring
        if (nnz(diff(pow)) == 0) & length(pow) == 3
            x = recover(pow(1));
            y = sdpvar(1,1);
            w = sdpvar(1,1);
            dx = double(x);
            assign([w y],[dx^3 dx^2]);
            DefinedMonoms(end+1).power = [pow(1:2)];
            DefinedMonoms(end).variable = getvariables(y);
            DefinedMonoms(end+1).power = [pow(3) getvariables(y)];
            DefinedMonoms(end).variable = getvariables(w);
            F = F + set(x*x == y) + set(x*y == w); 
            return
        end
        i = 1;
        found = 0;
        while ~found & i <= length(DefinedMonoms)
            if isequal(pow,DefinedMonoms(i).power)           
                w = recover(DefinedMonoms(i).variable);
                return
            end
            i = i + 1;
        end
        pow1 = pow(1:floor(length(pow)/2));
        pow2 = pow(1+floor(length(pow)/2):end);
        [w1,F,DefinedMonoms] = recdef(pow1,F,DefinedMonoms,setinitials);
        [w2,F,DefinedMonoms] = recdef(pow2,F,DefinedMonoms,setinitials);
        w = sdpvar(1,1);
        w1w2 = w1*w2;
        F = F + set(w1w2 == w);
        if setinitials
            assign(w,double(w1w2));
        end
        DefinedMonoms(end+1).power = pow;
        DefinedMonoms(end).variable = getvariables(w);
end