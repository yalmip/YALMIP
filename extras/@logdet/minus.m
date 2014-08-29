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

% function Z = minus(cx,P)
% %display           Overloaded
% 
% % Author Johan Löfberg 
% % $Id: minus.m,v 1.4 2007-02-07 09:11:27 joloef Exp $  
% 
% % Standard case c't-logdet(P)
% 
% if isa(P,'logdet') % sdpvr - logdet
% 
%     if prod(size(cx))>1
%         error('Only scalar terms can be added to a logdet term');
%     end
%     
%     if isa(cx,'logdet')
%         error('Logdet objects can only be added');
%     end
%     
%     Z = P;
%     if isempty(P.cx)
%         Z.cx = cx;
%     else
%         Z.cx = cx-P.cx;
%     end
%     Z.gain = -Z.gain;
% else % logdet - cx
%     temp = cx;
%     cx = P;
%     P = temp;
%     
%     Z = P;
%     if isempty(P.cx)
%         Z.cx = -cx;
%     else
%         Z.cx = P.cx-cx;
%     end
% end