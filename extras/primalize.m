function [Fdual,objdual,y,X] = primalize(F,obj)
% PRIMALIZE Create the dual of an SDP given in dual form
%
% [Fd,objd,y] = primalize(F,obj)
%
% Input
%  F   : Primal constraint in form C-Ay > 0, Fy = g
%  obj : Primal cost maximize b'y
%
% Output
%  Fd  : Dual constraints in form X>0, Trace(AiX)==bi+dt
%  obj : Dual cost trace(CX)
%  y   : The detected primal free variables
%
% Example
%  See the HTML help.
%
% See also DUAL, SOLVESDP, SDPVAR, DUALIZE

err = 0;

if isa(F,'constraint')
    F = (F);
end

% It's general, but not insanely general...
if ~(islinear(F) & islinear(obj))
    if nargout == 6
        Fdual = ([]);objdual = [];y = []; err = 1;
    else
        error('Can only primalize linear problems');
    end
end
if any(is(F,'socc'))
    if nargout == 6
        Fdual = ([]);objdual = [];y = []; err = 1;
    else
        error('Cannot primalize second order cone constraints');
    end
end
if isa(obj,'sdpvar')
    if any(is(F,'complex')) | is(obj,'complex')
    if nargout == 6
        Fdual = ([]);objdual = [];y = []; X = []; t = []; err = 1;
    else
        error('Cannot primalize complex-valued problems');
    end
    end
end
if any(is(F,'integer')) | any(is(F,'binary'))
    if nargout == 6
        Fdual = ([]);objdual = [];y = []; err = 1;
    else
        error('Cannot primalize discrete problems');
    end
end

% Create model using the standard code
[model,~,diagnostic] = export(F,obj,sdpsettings('solver','sedumi'),[],[],1);
if isempty(model)
    if isfield(diagnostic,'problem') && isfield(diagnostic,'info');
        error("export failed: model is empty, problem=%d,info=%s", ...
            diagnostic.problem,diagnostic.info);
    else
        error("export failed: model is empty");
    end
end

Fdual = ([]);
xvec = [];
if any(model.K.f)
    t = sdpvar(model.K.f,1);
    xvec = [xvec;t];
end

if any(model.K.l)
    x = sdpvar(model.K.l,1);
    xvec = [xvec;x];
    Fdual = Fdual + (x>=0);
end

if any(model.K.q)
    for i = 1:length(model.K.q)
        x = sdpvar(model.K.q(i),1);
        xvec = [xvec;x];
        Fdual = Fdual + (cone(x(2:end),x(1)));
    end
end

if any(model.K.s)
    for i = 1:length(model.K.s)
        X{i} = sdpvar(model.K.s(i),model.K.s(i));
        xvec = [xvec;X{i}(:)];
        Fdual = Fdual + (X{i}>=0);       
    end
end

objdual = model.C(:)'*xvec;
Fdual = Fdual + (-model.b == model.A'*xvec);

yvars = union(getvariables(F),getvariables(obj));
y = recover(yvars);

yalmip('associatedual',getlmiid(Fdual(end)),y);
