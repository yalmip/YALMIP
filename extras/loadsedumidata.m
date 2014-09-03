function [F,h] = loadsedumidata(varargin)
%LOADSEDUMIDATA Loads a problem definition in the SeDuMi format
%
%    [F,h] = loadsedumidata('filename')  Loads the problem min(h(x)), F(x)>0 from file 'filename'
%    [F,h] = loadsedumidata              An "Open" - box will be opened

filename = varargin{1};

% Does the file exist
if ~exist(filename)
    filename = [filename '.mat'];
    if ~exist(filename)
        error(['No such file.']);
    end
end
   
load(filename)
try   
    if ~exist('At')
        At = A;
    end
    if ~exist('b')
        b = zeros(size(At,1),1);
    else
        b = b(:);
    end
    if ~exist('c')
        if exist('C')
            c = C(:);
        else
            c = zeros(size(At,2),1);
        end
    else
        c = c(:);
    end
    
    K = K;
catch
    error('The file should contain the data At, b, c and K');
end

nvars = length(b);
x = sdpvar(nvars,1);

% No reason to try to do factor tracking here
x = flush(x);

if size(At,2)~=length(b)
    At = At';
end

F = ([]);
top = 1;

if isvalidfield(K,'f')
    X = c(top:top+K.f-1)-At(top:top+K.f-1,:)*x;
    F = F + (X(:) == 0);
    top = top + K.f;
end

if isvalidfield(K,'l')
    X = c(top:top+K.l-1)-At(top:top+K.l-1,:)*x;
    F = F + (X(:)>=0);
    top = top + K.l;
end

if isvalidfield(K,'q')
    for i = 1:length(K.q)
        X = c(top:top+K.q(i)-1)-At(top:top+K.q(i)-1,:)*x;
        F = F + (cone(X(2:end),X(1)));
        top = top + K.q(i);
    end
end

if isvalidfield(K,'r')
    for i = 1:length(K.r)
        X = c(top:top+K.r(i)-1)-At(top:top+K.r(i)-1,:)*x;
        F = F + (rcone(X(3:end),X(2),X(1)));
        top = top + K.r(i);
    end
end

if isvalidfield(K,'s')
    for i = 1:length(K.s)
        [ix,iy,iv] = find([c(top:top+K.s(i)^2-1) At(top:top+K.s(i)^2-1,:)]);
        off = (ix-1)/(K.s(i)+1);
        if all(off == round(off))
            X = c(top:top+K.s(i)^2-1)-At(top:top+K.s(i)^2-1,:)*x;
            if isa(X,'sdpvar')
                F = F + (diag(reshape(X,K.s(i),K.s(i))) >= 0);
            else
                X
                i
                'silly data!'
            end
            top = top + K.s(i)^2;
        else
            X = c(top:top+K.s(i)^2-1)-At(top:top+K.s(i)^2-1,:)*x;
            X = reshape(X,K.s(i),K.s(i));
            X = (X+X')/2;
            F = F + (X >= 0);
            top = top + K.s(i)^2;
        end
    end
end

h = -b'*x;

function ok = isvalidfield(K,fld)
ok = 0;
if isfield(K,fld)
    s = getfield(K,fld);
    if prod(size(s))>0
        if s(1)>0
            ok = 1;
        end
    end
end

