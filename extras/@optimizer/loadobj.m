function P = loadobj(P)
if isstruct(P),
    P = class(P,'optimizer');
end