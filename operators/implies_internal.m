function varargout = implies_internal(varargin)

X = varargin{1};
Y = varargin{2};

if nargin == 2
    zero_tolerance = 1e-5;
else
    zero_tolerance = varargin{3};
end

% Normalize
if isa(X,'constraint')
    X = lmi(X,[],[],1);
end
if isa(Y,'constraint')
    Y = lmi(Y,[],[],1);
end

% % Special case x == 1 (or x = 0) implies something
% Avoids introducing a glue binary for the equality
% Fixes #690
if isa(X,'lmi') && length(X)==1 && is(X,'equality')
    x = sdpvar(X);
    if is(x,'binary')
        B = getbase(x);
        if isequal(B,[1 -1]) || isequal(B,[-1 1])
            varargout{1} = implies_internal(recover(depends(x)),Y);
            return
        elseif isequal(B,[0 -1]) || isequal(B,[0 1])
            varargout{1} = implies_internal(1-recover(depends(x)),Y);
            return
        end        
    end
end

if isa(X,'sdpvar') & isa(Y,'sdpvar')
    varargout{1} = (Y >= X);
elseif isa(X,'sdpvar') & isa(Y,'lmi')
    varargout{1} = binary_implies_constraint(X,Y);
elseif isa(X,'lmi') & isa(Y,'sdpvar')
    varargout{1} = constraint_implies_binary(X,Y,zero_tolerance);
elseif isa(X,'lmi') & isa(Y,'lmi')
    % Glue them using a binary
    binvar d
    varargout{1} = [constraint_implies_binary(X,d,zero_tolerance), binary_implies_constraint(d,Y)];
end

function F = binary_implies_constraint(X,Y)
if isempty(Y) | length(Y)==0
    F = [];
    return
end
switch settype(Y)
    case 'multiple'
        % recursive call if we have a mixture
        Y1 = Y(find(is(Y,'equality')));
        Y2 = Y(find(is(Y,'elementwise')));
        Y3 = Y(find(is(Y,'socp')));
        Y4 = Y(find(is(Y,'sdp')));
        F = [binary_implies_constraint(X,Y1),
            binary_implies_constraint(X,Y2),
            binary_implies_constraint(X,Y3),
            binary_implies_constraint(X,Y4)];
        
    case 'elementwise' % X --> (Y(:)>=0)
        Y = sdpvar(Y);
        Y = Y(:);
        F = binary_implies_linearnegativeconstraint(-Y,X);
        
    case 'equality'    % X --> (Y(:)==0)
        Y = sdpvar(Y);
        Y = Y(:);
        [M,m,infbound]=derivebounds(Y);
        if infbound
            warning('You have unbounded variables in IMPLIES leading to a lousy big-M relaxation. Read more https://yalmip.github.io/faq/warningbigm/');                        
        end
        F = binary_implies_linearequality(Y,reshape(X,[],1),M,m);
        % F = [F, binary_implies_linearnegativeconstraint(-Y,reshape(X,[],1),-m,-M)];
        
    case 'sdp'         % X --> (Y>=0)
        if length(X)>1
            error('IMPLIES not implemented for vector x implies lmi.');
        end
        Y = sdpvar(Y);
        % Elements in matrix
        y = Y(find(triu(ones(length(Y)))));
        % Derive bounds on all elements
        [M,m,infbound]=derivebounds(y);
        if infbound
            warning('You have unbounded variables in IMPLIES leading to a lousy big-M relaxation. Read more https://yalmip.github.io/faq/warningbigm/');            
        end
        % Crude lower bound eig(Y) > -max(abs(Y(:))*n*I
        m=-max(abs([M;m]))*length(Y);
        % Big-M relaxation...
        F = [Y >= (1-X)*m*eye(length(Y))];
        
    otherwise
        error('IMPLIES only implemented for linear (in)equalities and semidefinite constraints');
end

function F = constraint_implies_binary(X,Y,zero_tolerance)

switch settype(X)
    case 'multiple'
        X1 = X(find(is(X,'equality')));
        X2 = X(find(is(X,'elementwise')));
        if length(X1)+length(X2) < length(X)
            disp('Binary variables can only be activated by linear (in)equalities in IMPLIES');
        end
        d1 = binvar(length(Y),1);
        d2 = binvar(length(Y),1);
        F1 = constraint_implies_binary(X1,d1,zero_tolerance);
        F2 = constraint_implies_binary(X2,d2,zero_tolerance);
        F = [F1, F2, Y >= d1+d2-1];
        
    case 'elementwise'
        X = -sdpvar(X);
        
        if length(Y)==length(X)
            % Elementwise implies, either by user or internally
            F = linearnegativeconstraint_implies_binary(X,Y,[],[],zero_tolerance);
        elseif length(Y)== 1
            % Many rows should imply one. Create intermediate and use AND
            d = binvar(length(X),1);
            F = [linearnegativeconstraint_implies_binary(X,d,[],[],zero_tolerance), Y >= sum(d)+1-length(X)];
        else
            error('Inconsistent sizes in implies_internal')
        end
        
    case 'equality'
        X = sdpvar(X);X = reshape(X,[],1);
        n = length(X);
        if 0 % isequal(getbase(X),[-ones(n,1) eye(n)]) | isequal(getbase(X),[ones(n,1) -eye(n)]) & all(ismember(depends(X),yalmip('binvariables')));
            % Smart code for X == 1 implies Y
            F = [Y(:) >= recover(getvariables(X))];
        else           
            d = binvar(length(X),2,'full');
            %F = linearnegativeconstraint_implies_binary([-zero_tolerance-X;X-zero_tolerance],[d(:,1);d(:,2)],[],[],zero_tolerance);
            F = linearnegativeconstraint_implies_binary([-zero_tolerance-X;X-zero_tolerance],[d(:,1);d(:,2)],[],[],zero_tolerance/100);
            if length(X)==length(Y)
                % elementwise version
                F = [F, Y >= d(:,1)+d(:,2)-1];
            elseif length(Y)==1
                % All elements must be zero
                F = [F, Y >= sum(sum(d))-2*length(X)+1];
            else
                error('Inconsistent sizes in implies_internal')
            end
        end
        
    otherwise
        disp('IMPLIES not implemented for this case');
        error('IMPLIES not implemented for this case');
end



