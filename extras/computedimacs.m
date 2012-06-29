function dimacs = computedimacs(b,c,A,xin,y,s,K);
% COMPUTEDIMACS
%
% min <C,X> s.t     AX = b, X > 0
% max b'y   s.t S-C+A'y =0, S > 0

% If no primal exist, fake till later
if isempty(xin)
    x = c*0;
else
    x = xin;
end

if isempty(s)
    s = c-A'*y;
end

xres = inf;
sres = inf;

% Not officially defined in DIMACS
if K.f>0
    sres = -min(norm(s(1:K.f),inf));
end

% Errors in linear cone
if K.l>0
    xres = min(x(1+K.f:K.f+K.l));
    sres = min(s(1+K.f:K.f+K.l));
end

% Errors in quadratic cone
if K.q(1)>0
    top = K.f+K.l;
    for i = 1:length(K.q)
        X = x(1+top:top+K.q(i));
        S = s(1+top:top+K.q(i));
        xres = min(xres,X(1)-norm(X(2:end)));
        sres = min(sres,S(1)-norm(S(2:end)));
        top = top + K.q(i);
    end
end

% Errors in semidefinite cone
if K.s(1)>0
    top = K.f+K.l+K.q+K.r;
    for i = 1:length(K.s)
        X = reshape(x(1+top:top+K.s(i)^2),K.s(i),K.s(i));
        S = reshape(s(1+top:top+K.s(i)^2),K.s(i),K.s(i));
        xres = min(xres,min(eig(full(X))));
        sres = min(sres,min(eig(full(S))));
        top = top + K.s(i)^2;
    end
end

err1 = norm(b-A*x)/(1 + norm(b,inf));
err2 = max(0,-xres)/(1 + norm(b,inf));
err3 = conenorm(s-(c-A'*y),K)/(1+norm(c,inf));
%err3 = norm(s-(c-A'*y))/(1+norm(c,inf)); % Used by some solvers
err4 = max(0,-sres)/(1+max(abs(c)));
err5 = (c'*x-b'*y)/(1+abs(c'*x)+abs(b'*y));
err6 = x'*(c-A'*y)/(1+abs(c'*x)+abs(b'*y));

% No primal was computed
if isempty(xin)
    err1 = nan;
    err2 = nan;
    err5 = nan;
    err6 = nan;
end
dimacs = [err1 err2 err3 err4 err5 err6];

function t = conenorm(s,K)

% Implementation of the norm described on
% http://plato.asu.edu/dimacs/node3.html

t = 0;

if K.f + K.l>0
    t = t + norm(s(1:K.f+K.l));
end

top = 1+K.f+K.l;
if K.q(1)>0
    for i = 1:length(K.q)
        t = t + norm(s(top:top+K.q(i)-1));
        top  = top + K.q(i);
    end
end

if K.s(1)>0
    for i = 1:length(K.s)
        S = reshape(s(top:top+K.s(i)^2-1),K.s(i),K.s(i));
        t = t + norm(S,'fro');
        top  = top + K.s(i)^2;
    end
end