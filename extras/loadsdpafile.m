function [F,h] = loadsdpafile(varargin)
%LOADSDPAFILE Loads a problem definition in the SDPA format
%
%    [F,h] = loadsdpafile('filename')  Loads the problem min(h(x)), F(x)>0 from file 'filename'
%    [F,h] = loadsdpafile         A "Open" - box will be opened

filename = varargin{1};

% Does the file exist
if ~exist(filename)
    filename = [filename '.dat-s'];
    if ~exist(filename)
        error(['No such file.']);
    end
end
    
% Load using SeDuMi
try
    [At,b,c,K] = fromsdpa(filename);
catch
    error('LOADSDPAFILE currently requires SeDuMi to be installed');
end

nvars = length(b);
x = sdpvar(nvars,1);

F = ([]);
top = 1;
if isfield(K,'l')
    if K.l(1)>0
        X = c(top:top+K.l-1)-At(top:top+K.l-1,:)*x;
        F = F + (X(:)>=0);
        top = top + K.l;
    end
end

if isfield(K,'s')
    if K.s(1)>0
        for i = 1:length(K.s)
            [ix,iy,iv] = find([c(top:top+K.s(i)^2-1) At(top:top+K.s(i)^2-1,:)]);
            off = (ix-1)/(K.s(i)+1);
            if all(off == round(off))           
                X = c(top:top+K.s(i)^2-1)-At(top:top+K.s(i)^2-1,:)*x;
                F = F + (diag(reshape(X,K.s(i),K.s(i))) >= 0);
                top = top + K.s(i)^2;
            else
                X = c(top:top+K.s(i)^2-1)-At(top:top+K.s(i)^2-1,:)*x;
                F = F + (reshape(X,K.s(i),K.s(i)) >= 0);
                top = top + K.s(i)^2;
            end
        end
    end
end

h = -b'*x;

