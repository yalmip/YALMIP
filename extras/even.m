function YESNO = even(x,order,d)
%EVEN Check if all numbers are even

% Author Johan Löfberg
% $Id: even.m,v 1.2 2004-07-02 08:17:30 johanl Exp $

if (nargin==2) | (nargin==3)
    switch order 
        case 'rows'
           YESNO = ~any(rem(x,d),2);         
        case 'cols'
            YESNO = ~any(rem(x,d),1);
        otherwise
            error
    end
else
YESNO = nnz(rem(x,2))==0;
end

