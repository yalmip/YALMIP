function Z = plus(A,B)
%plus           Overloaded

% Author Johan Löfberg 
% $Id: plus.m,v 1.5 2007-02-07 09:11:27 joloef Exp $  

% Standard case A.cx + A.logdetP + (B.cx + B.logdetP)

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

if isequal(Acx+Bcx,0)
    Z.cx = [];
else
    Z.cx = Acx + Bcx;
end
    
Z.P = {Alog{:},Blog{:}};
Z.gain = [Again Bgain];

% function Z = plus(X,Y)
% %display           Overloaded
% 
% % Author Johan Löfberg 
% % $Id: plus.m,v 1.5 2007-02-07 09:11:27 joloef Exp $  
% 
% % LOGDET + SDPVAR
% if isa(Y,'logdet')
%     Z = X;
%     X = Y;
%     Y = Z;
% end
% 
% if prod(size(Y))>1
%     error('Only scalar terms can be added to a logdet term');
% end
% 
% if isa(Y,'logdet')
%     Z = X;
%     Z.P = {Z.P{:},Y.P{:}};%blkdiag(Z.P,Y.P);
%     Z.gain = [X.gain Y.gain];
%     return
% end
% 
% Z = X;
% if isempty(Z.cx)
%     Z.cx = Y;
% else
%     Z.cx = plus(Z.cx,Y);
% end