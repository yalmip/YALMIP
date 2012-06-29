function d = det(X)
%DET (overloaded)

% Author Johan Löfberg 
% $Id: det.m,v 1.1 2006-08-10 18:00:19 joloef Exp $  

n = X.dim(1);
m = X.dim(2);

if n~=m
    error('Matrix must be square.');
else
    switch n
        case 1
            d = X;
        case 2
            % Freakin overloading on multiplication doesn't work. Probalby
            % stupid code...
            Y1.type = '()';
            Y2.type = '()';
            Y3.type = '()';
            Y4.type = '()';
            Y1.subs = {1,1};
            Y2.subs = {2,2};
            Y3.subs = {1,2};
            Y4.subs = {2,1};
            d = subsref(X,Y1)*subsref(X,Y2)-subsref(X,Y3)*subsref(X,Y4);
        otherwise
            d = 0;
            Y.type = '()';
            for i = 1:n   
                Y.subs = {i,1};
                xi = subsref(X,Y);
                if ~isequal(xi,0)
                    Y.subs = {[1:1:i-1 i+1:1:n],2:n};
                    subX = subsref(X,Y);
                    d = d + (-1)^(i+1)*xi*det(subX);
                end
            end
    end
end
% Reset info about conic terms
d.conicinfo = [0 0];