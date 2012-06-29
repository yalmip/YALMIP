function Z = plus(A,B)
%plus           Overloaded

% Author Johan Löfberg 
% $Id: mtimes.m,v 1.2 2007-10-03 12:31:41 joloef Exp $  

% Standard case A.cx + A.logdetP + (B.cx + B.logdetP)

if ~(isa(A,'double') | isa(B,'double'))
    error('LOGDET objects can only be multiplied with constants')
end

if isa(A,'logdet')
    temp = A;
    A = B;
    B = temp;
end

% OK, constant*logdet
% end
B.gain = B.gain*A;
B.cx = B.cx*A;
Z = B;