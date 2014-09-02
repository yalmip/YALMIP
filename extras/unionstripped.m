function c = unionstripped(a,b)
%UNIONSTRIPPED  Internal function (version without checkings etc.)

c = uniquestripped([a(:)' b(:)']);


