function [Fimage,objimage,x,y] = imagemodel(F,obj,options)
% IMAGEMODEL Explicitely removes equality constraints from model.
%
% [Fi,hi,x,y] = imagemodel(F,h)
%
% Input
%  F   : Constraint in form F(x)>0, Ax=b
%  h   : objective function h(x)
%
% Output
%  Fi  : Constraints in form F(y)>0
%  hi  : objective function h(y)
%  x   : original variables
%  z   : old variables expressed in basis y
%
% To obtain a solution in the original variables x, use assign(x,double(y))
%
% Note that this reduction is automatically done when you call solvesdp and use 
% a solver that cannot handle equality constraints. Hence, there is typically no
% reason to use this command, unless some further manipulations are going
% to be done.
%
% See also DUALIZE, PRIMALIZE

% Check for unsupported problems
err = 0;
p1 = ~(isreal(F) & isreal(obj));
p2 = ~(islinear(F) & islinear(obj));
p3 = any(is(F,'integer')) | any(is(F,'binary'));
if p1 | p2 | p3
    if nargout == 5
        Fdual = ([]);objdual = [];y = []; X = []; t = []; err = 1;
    else
        problems = {'Cannot imagalize complex-valued problems','Cannot imagalize nonlinear problems','Cannot imagalize discrete problems'};
        error(problems{min(find([p1 p2 p3]))});
    end
end

if nargin < 3
    options = sdpsettings('solver','sedumi','remove',1);
end
    
if any(is(F,'equality'))
    [model,recoverdata,solver,diagnostic,F] = compileinterfacedata(F,[],[],obj,options,0);
    if isfield(diagnostic,'problem')
        if diagnostic.problem == 1
            warning('Problem is infeasible');
            Fimage = [];
            objimage = [];
            x = [];
            y = [];
            return
        end
    end
    if isempty(model)
            Fimage = [];
            objimage = [];
            x = [];
            y = [];
            warning('Reduced problem does not have free variables. Optimal solution computed');
            return
    end
else
    Fimage = F;
    objimage = obj;
    x  = recover(unique([depends(obj) depends(F)]));
    y = x;
    return
end

y = sdpvar(length(model.c),1);

vecF = model.F_struc*[1;y];
K = model.K;
Fimage = ([]);
start = 1;
if model.K.l > 0
    Fimage = Fimage + (vecF(start:start+K.l-1) >= 0);
    start = start + K.l;
end

if model.K.q(1) > 0
    for i = 1:length(model.K.q)
        z = vecF(start:start+K.q(i)-1)
        Fimage = Fimage + (cone(z(2:end),z(1)));
        start = start + K.q(i);
    end
end

if model.K.s(1)>0
    for i = 1:length(model.K.s)
        z = vecF(start:start+K.s(i)^2-1);
        Fimage = Fimage + (reshape(z,K.s(i),K.s(i)));
        start = start + K.s(i)^2;
    end
end

objimage = model.c'*y+y'*model.Q*y+model.f;

x = recover(recoverdata.used_variables);
y = recoverdata.H*y+recoverdata.x_equ;
