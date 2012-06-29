function YESNONA = isconvex(p)
%ISCONVEX Checks if scalar function is convex

% Author Johan Löfberg 
% $Id: isconvex.m,v 1.2 2005-10-05 20:50:42 joloef Exp $   
p=p;
if is(p,'linear')
    YESNONA = 1;
    return;
elseif is(p,'quadratic')
    [Q,c,f,x,info] = quaddecomp(p);
    if ~info
        if all(real(eig(Q+Q')) > -1e-13)
            YESNONA = 1;
        end
    end
end

vars = depends(p);
x = recover(depends(p));
convex = 1;
iterations = 0;
while convex & iterations<10
    y1 = randn(length(vars),1);
    assign(x,y1);
    p1 = double(p);
    y2 = randn(length(vars),1);
    assign(x,y2);
    p2 = double(p);
    yc = ((y1+y2)/2);
    assign(x,yc);
    pc = double(p);
    if pc>(p1+p2)/2
        convex = 0;
    end
    iterations = iterations + 1;
end

% Maybe we didn't manage to prove non-convexity
if convex
    H = hessian(p,x);
    v = sdpvar(length(H),1);
    sol = solvesos(set(sos(v'*H*v)),[],sdpsettings('verbose',1));
    if sol.problem == 0
        YESNONA = 1;
    else
        YESNONA = nan;
    end       
else
    YESNONA = 0;
end

