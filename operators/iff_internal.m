function varargout = iff_internal(varargin)
X = varargin{1};
Y = varargin{2};

if nargin == 2
    zero_tolerance = 1e-5;
else
    zero_tolerance = abs(varargin{3});
end

% Normalize data
if isa(Y,'constraint')
    Y=lmi(Y,[],[],1);
end
if isa(X,'constraint')
    X=lmi(X,[],[],1);
end
if isa(X,'lmi') & isa(Y,'sdpvar')
    temp = X;
    X = Y;
    Y = temp;
end

if isa(X,'sdpvar') & isa(Y,'sdpvar')
    % user presumably works with binary variables
    varargout{1} = [Y == X];
    
elseif (isa(X,'sdpvar') & isa(Y,'lmi'))
    % Binary variable iff a constraint holds
    switch settype(Y)
        case 'elementwise'  % binary X <--> Y(:)>=0
            varargout{1} = binary_iff_lp(X,-sdpvar(Y),zero_tolerance);
        case 'equality'
            varargout{1} = binary_iff_equality(X,sdpvar(Y),zero_tolerance);
        case 'multiple'
            Y1 = Y(find(is(Y,'equality')));
            Y2 = Y(find(is(Y,'elementwise')));
            % Glue using intermediate variables holding truth of each set
            % of constraints. X is true iff both di are true and v.v.
            % Hence, we make a recursive call with two new models
            binvar d1 d2
            C1 = iff_internal(d1,Y1,zero_tolerance);
            C2 = iff_internal(d2,Y2,zero_tolerance);
            C3 = [X <= d1,X <= d2,X >= d1+d2-1];
            varargout{1} = [C1, C2, C3];
          %  varargout{1} = [iff_internal(X,Y1),iff_internal(X,Y2)];
        otherwise
            error('IFF not implemented for this case');
    end    
elseif isa(X,'lmi') & isa(Y,'lmi')
    % Constraint holds iff constraint holds. 
    % Glue using an intermediate variable
    binvar d    
    C1 =  iff_internal(d,X,zero_tolerance);
    C2 =  iff_internal(d,Y,zero_tolerance);
    varargout{1} = [C1,C2];
end

function F = binary_iff_lp(X,f,zero_tolerance)
% X == 1    <=>   f<=0
[M,m,infbound] = derivebounds(f);
if infbound
    warning('You have unbounded variables in IFF leading to a lousy big-M relaxation.');
end

[nf,mf]=size(f);
if nf*mf==1    
    F = linearnegativeconstraint_iff_binary(f,X,M,m,zero_tolerance);
else
    f = reshape(f,nf*mf,1);
    di = binvar(nf*mf,1);
    F = linearnegativeconstraint_iff_binary(f,di,M,m,zero_tolerance);
    if length(X)==1
        % X is true if any di
        F = [F, X>=sum(di)-length(di)+1, X <= di];
    else
        % This must be a vectorized X(i) iff f(i)
        F = [F, di == X];
    end
    
    % di=0 means the ith hypeplane is violated
    % X=1 means we are in the polytope
    % F  = [f <= M*(1-X), f>=eps+(m-eps).*di, X>=sum(di)-length(di)+1, X <= di];

    % Add some cuts for c < a'x+b < d
    [bA] = getbase(f);
    b = bA(:,1);
    A = bA(:,2:end);
    S = zeros(0,length(di));
    for i = 1:length(b)
        j = findrows(abs(A),abs(A(i,:)));
        j = j(j > i);
        if length(j)==1
            S(end+1,[i j]) = 1;
        end
    end
    if size(S,1) > 0
        % Add cut cannot be outside both constraints
        F = F + (S*di >= 1);
    end
end

function F = binary_iff_equality(X,Y,zero_tolerance)

% Things like iff(x,y==1) is sent as X, 1-Y
if isLogicalVector(X) && isLogicalVector(Y)
    if isequal(getbase(Y),[0 -1]) | isequal(getbase(Y),[-1 1])
        Y = -Y;
    end
     F = [1-X == Y];
     return
 end
Y = Y(:);
d = binvar(length(Y),3);
% We have to model every single line of equality.
% di1 : Yi < -eps
% di2 : -eps < Yi < eps
% di3 : eps < Yi
C = sum(d,2)==1;
[M,m,infbounds] = derivebounds(Y);
for i = 1:length(Y)
    if all(ismember(depends(Y(i)),yalmip('binvariables'))) & all(getbase(Y(i))==fix(getbase(Y(i))))
        eps = 0.5;
    else
        eps = 0;
    end
    C1 =  binary_iff_lp(d(i,1),Y(i)+eps,abs(zero_tolerance));             % Y <-eps
    C2 =  binary_iff_lp(d(i,2),[Y(i)-eps;-eps-Y(i)],abs(zero_tolerance)); % -eps < Y < eps
    C3 =  binary_iff_lp(d(i,3),eps-Y(i),abs(zero_tolerance));             % eps < Y 
    C = [C, C1,C2,C3];
    % Cut off some trivial cases
    if m(i)>=0
        C = [C,d(i,1)==0];
    elseif M(i)<=0
        C = [C,d(i,3)==0];
    end
end
% X is true if all di2 are true and v.v
F = [C,X <= d(:,2), X >= sum(d(:,2))-length(Y)+1];

function yes = isLogicalVector(X)
yes = 0;
X = sdpvar(X);
B = getbase(X);
if all(ismember(depends(X),yalmip('binvariables')))
    b = B(:,1);
    if all(b==0 | b==1)
        A = B(:,2:end);
        if all(A == 0 | A == 1 | A == -1)
            if all(sum(A | A,2)==1)
                if all(sum(B,2)<= 1)
                    yes = 1;
                end
            end
        end
    end        
end