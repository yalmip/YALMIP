function y = power(x,d)
%POWER (overloaded)

% Author Johan Löfberg 
% $Id: power.m,v 1.13 2009-10-14 07:28:42 joloef Exp $   


% Vectorize x if d is vector
if prod(size(x))==1 & (prod(size(d))>1)
    x = x.*ones(size(d));
end

x = flush(x);

% Reuse code
if prod(size(x))==1 & (prod(size(d))==1)
    y = mpower(x,d);
    return 
end

% Sanity Check
if prod(size(d))>1
    if any(size(d)~=size(x))
        error('Matrix dimensions must agree.');
    end
else
    d = ones(x.dim(1),x.dim(2))*d;
end

% Trivial cases
if isa(d,'double')
    if all(all(d==0))
        if x.dim(1)~=x.dim(2)
            error('Matrix must be square.')
        end
        y = eye(x.dim(1),x.dim(2)).^0;
        return
    end
    if all(all(d==1))
        y = x;
        return
    end
end

% Fractional, negative or different powers are
% treated less efficiently using simple code.
fractional = any(any((ceil(d)-d>0)));
negative = any(any(d<0));
different = ~all(all(d==d(1)));
if fractional | negative | different
    if x.dim(1)>1 | x.dim(2)>1
        [n,m] = size(x);        
        y = [];
        for i = 1:n % FIX : Vectorize!
            temp = [];
            for j = 1:m
                temp = [temp extsubsref(x,i,j).^d(i,j)];
            end
            y = [y;temp];
        end
        return
    else
        base = getbase(x);
        if isequal(base,[0 1])
            mt = yalmip('monomtable');
            var = getvariables(x);
            previous_var = find((mt(:,var)==d)  & (sum(mt~=0,2)==1));
            if isempty(previous_var)
                mt(end+1,:) = mt(getvariables(x),:)*d;
                yalmip('setmonomtable',mt);
                y = recover(size(mt,1));
            else
                y = recover(previous_var);
            end
        elseif (size(base,2) == 2) & base(1)==0
            % Something like a*t^-d
            y = base(2)^d*recover(getvariables(x))^d;
        else
            error('Only unit scalars can have negative or non-integer powers.');
        end
    end
    return
end

% Back to scalar power...
d = d(1,1);
if x.dim(1)>1 | x.dim(2)>1
    switch d
        case 0
            y = 1;
        case 1
            y = x;
        otherwise
            y = x.*power(x,d-1);
    end
else
    error('This should not appen. Report bug (power does not use mpower)')
end