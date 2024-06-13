function Z = minus(A,B)
%minus           Overloaded

% Standard case A.cx + A.logdetP - (B.cx + B.logdetP)

if isa(A,'sdpvar') | isa(A,'double') 
    if prod(size(A))>1
        error('Only scalar terms can be added to a logdet term');
    end
    Acx = A;
    Alog = {};
    Again = [];
else
    Z = A;
    Acx  = A.cx;
    Alog = A.P;
    Again = A.gain;
end

if isa(B,'sdpvar') | isa(B,'double') 
    if prod(size(B))>1
        error('Only scalar terms can be added to a logdet term');
    end
    Bcx = B;
    Blog = {};
    Bgain = [];
else
    Z = B;
    Bcx  = B.cx;
    Blog = B.P;
    Bgain = B.gain;
end

if isempty(Acx)
    Acx = 0;
end
if isempty(Bcx)
    Bcx = 0;
end

if isequal(Acx-Bcx,0)
    Z.cx = [];
else
    Z.cx = Acx - Bcx;
end
    
Z.P = {Alog{:},Blog{:}};
Z.gain = [Again -Bgain];

