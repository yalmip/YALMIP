function YESNO = even(x,order,d)
%EVEN Check if all numbers are even

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

