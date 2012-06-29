function p = ellipsoid(X)
% polytope  Converts set object to ELLIPSOID object        
%
% E     = ellipsoid(F)
% [P,x] = ellipsoid(F)
%
% P : ELLIPSOID object (Requires the Ellipsoidal Toolbox)
% x : sdpvar object defining the variables in the polytope P.H*x<P.K

% Author Johan Löfberg
% $Id: ellipsoid.m,v 1.1 2005-10-02 19:12:14 joloef Exp $


if all(is(X,'element-wise'))% & all(is(X,'linear'))
    f = [];
    for i = 1:length(X)
        if  X.clauses{i}.type==2
            fi =  X.clauses{i}.data;
            [Q,c,f,x,info] = quaddecomp(fi);       
            p=ellipsoid(c,-Q);
        end
    end
    
else
    error('ELLIPSOID can only be applied to SET objects with quadratic inequalities.')
end