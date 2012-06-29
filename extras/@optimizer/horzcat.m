function P = horzcat(A,B)

% Author Johan Löfberg
% $Id: horzcat.m,v 1.11 2007-08-17 13:17:16 joloef Exp $

if isa(A,'optimizer') & isa(B,'optimizer')
    if isequal(A.input.expression,B.input.expression) & isequal(A.output.expression,B.output.expression)
        P = optimizer([A.F,B.F],A.h+B.h,A.ops,A.input.expression,A.output.expression);
    else
        error('The optimizer objects must share parameters and decision variables');
    end
elseif isa(B,'optimizer') & (isa(A,'constraint') | isa(A,'lmi'))
    P = A;
    A = B;
    B = P;
end

P = optimizer([A.F,B],A.h,A.ops,A.input.expression,A.output.expression);


