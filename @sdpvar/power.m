function y = power(x,d)
%POWER (overloaded)

% Vectorize x if d is vector
if numel(x)==1 & (numel(d)>1)
    x = x.*ones(size(d));
end
% Vectorize if x is a vector
if numel(d)==1 & (numel(x)>1)
    d = d.*ones(size(x));
end
if ~isequal(size(d),size(x))
    error('Dimension mismatch in power');
end

% Reuse code
if numel(x)==1 && numel(d)==1
    y = mpower(x,d);
    return 
end

if isa(d,'sdpvar')
    % Call helper which vectorizes the elements
    y = powerinternalhelper(d,x);
    if isa(y,'sdpvar')
        y.extra.createTime = definecreationtime;
    end
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
if isnumeric(d)
    if all(all(d==0))
        y = ones(x.dim(1),x.dim(2));
        return
    end
    if all(all(d==1))
        y = x;
        return
    end
    if isnan(d)
        disp('You have NaNs in model (<a href="yalmip.github.io/naninmodel">learn to debug</a>)')
        error('NaN power makes no sense.');
    end
end

% Fractional, negative or different powers are
% treated less efficiently using simple code.
fractional = any(any((ceil(d)-d>0)));
negative = any(any(d<0));
different = ~all(all(d==d(1)));
if fractional | negative | different
    if x.dim(1)>1 | x.dim(2)>1
        if isequal(x.basis,[spalloc(prod(x.dim),1,0) speye(prod(x.dim))]) & all(d==d(1))
            % Simple case x.^d
            y = vectorizedUnitPower(x,d);
            y.extra.createTime = definecreationtime;
            return
        end
        [n,m] = size(x);        
        y = [];
        for i = 1:n % FIX : Vectorize!
            if m == 1
                temp = extsubsref(x,i,1).^d(i,1);
            else
                temp = [];
                for j = 1:m
                    temp = [temp extsubsref(x,i,j).^d(i,j)];
                end
            end
            y = [y;temp];
        end
        y.extra.createTime = definecreationtime;
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
    y.extra.createTime = definecreationtime;
    return
end

if isequal(x.basis,[spalloc(prod(x.dim),1,0) speye(prod(x.dim))]) & all(d==d(1))
     % Simple case x.^d
     y = vectorizedUnitPower(x,d);
     y.extra.createTime = definecreationtime;
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

function y = vectorizedUnitPower(x,d)
d = d(1);
[mt,variabletype,hashM,hash] = yalmip('monomtable');
var = getvariables(x);

usedmt = mt(var,:);
newmt = usedmt*d;
hashV = newmt*hash;
if ~any(ismember(hashV,hashM))      
    variabletype = [variabletype newvariabletypegen(newmt)];
    y = size(mt,1) + (1:length(var));
    mt = [mt;newmt];
else
    y = [];
    allnewmt = [];
    newvariables = 0;
    keep = zeros(size(usedmt,1),1);
    for i = 1:length(hashV)
        previous_var = find(abs(hashM - hashV(i)) < 1e-20);
        if isempty(previous_var)
            newmt =  usedmt(i,:)*d;           
            variabletype = [variabletype newvariabletypegen(newmt)];
            newvariables = newvariables + 1;
            keep(i) = 1;
            y = [y size(mt,1)+newvariables];
        else
            y = [y previous_var];
        end
    end
    mt = [mt;usedmt(find(keep),:)*d];
end
yalmip('setmonomtable',mt,variabletype);
y = reshape(recover(y),x.dim);

            
