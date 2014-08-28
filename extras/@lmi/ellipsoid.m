function p = ellipsoid(X)
% polytope  Converts set object to ELLIPSOID object        
%
% E     = ellipsoid(F)
% [P,x] = ellipsoid(F)
%
% P : ELLIPSOID object (Requires the Ellipsoidal Toolbox)
% x : sdpvar object defining the variables in the polytope P.H*x<P.K

X = flatten(F);
if all(is(X,'element-wise'))
    f = [];
    for i = 1:length(X)
        if  X.clauses{i}.type==2
            fi =  X.clauses{i}.data;
            [Q,c,f,x,info] = quaddecomp(fi);       
            p=ellipsoid(c,-Q);
        end
    end
    
else
    error('ELLIPSOID can only be applied to constraint objects with quadratic inequalities.')
end