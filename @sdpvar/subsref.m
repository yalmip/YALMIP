function varargout = subsref(varargin)
%SUBSREF (overloaded)

% Stupid first slice call (supported by MATLAB)
% x = sdpvar(2);x(1,:,:)
Y = varargin{2};
if length(Y)==1
    if  length(Y.subs) > 2 && isequal(Y.type,'()') && ~( isa(Y(1).subs{1},'constraint') || isa(Y(1).subs{1},'lmi') )
        i = 3;
        ok = 1;
        while ok && (i <= length(Y.subs))
            ok = ok && (isequal(Y.subs{i},1) || isequal(Y.subs{i},':'));
            i = i + 1;
        end
        if ok
            Y.subs = {Y.subs{1:2}};
        else
            error('??? Index exceeds matrix dimensions.');
        end
    end
end


X = varargin{1};

try
    switch Y(1).type
        case '()'
            if  isa(Y(1).subs{1},'constraint') || isa(Y(1).subs{1},'lmi') 
                z = Y(1).subs{1};
                for i = 1:length(Y(1).subs)
                    z = [z, Y(1).subs{i}];
                end
                % z == 0 => b+Ax == 0 => x = -A\b
                z = sdpvar(z);
                B = getbase(z);
                A = B(:,2:end);
                b = B(:,1);
                x0 = -A\b;
                x = recover(getvariables(z));
                varargout{1} = replace(X,x,x0);
                return                
            end
            % Check for simple cases to speed things up (yes, ugly but we all want speed don't we!)
            switch size(Y(1).subs,2)
                case 1
                    y = subsref1d(X,Y(1).subs{1},Y);                    
                case 2
                    y = subsref2d(X,Y.subs{1},Y(1).subs{2},Y);
                otherwise
                    if all( [Y(1).subs{3:end}]==1)
                        y = subsref2d(X,Y.subs{1},Y(1).subs{2},Y);
                    else
                        error('Indexation error.');
                    end
            end
        case '{}'
            varargout{nargout} = [];
            
            % it could be the case that we have an extended variable
            % This is a bit tricky, so we do the best we can; assume that
            % we want to replace the internal argument wih the new
            % expression
            OldArgument = recover(depends(X));
            vars = getvariables(X);
            mpt_solution = 1;
            if all(ismembc(vars,yalmip('extvariables')))
                for i = 1:length(X)
                    nonlinearModel = yalmip('extstruct',vars);
                    if isequal(nonlinearModel{1}.fcn,'pwa_yalmip') | isequal(nonlinearModel{1}.fcn,'pwq_yalmip')
                    else
                        mpt_solution = 0;
                    end
                end
                if mpt_solution
                    assign(nonlinearModel{1}.arg{2},Y(1).subs{:});
                    XX = value(X);
                    varargout{1} = value(X);
                    return
                end
            end
            vars = getvariables(X);
            if (length(vars) == 1) & ismembc(vars,yalmip('extvariables'))
                nonlinearModel = yalmip('extstruct',vars);
                OldArgument = [];
                for i = 1:length(nonlinearModel.arg)
                    if isa(nonlinearModel.arg{i},'sdpvar')
                        OldArgument = [OldArgument;  nonlinearModel.arg{i}];
                    end
                end
                if isnumeric([Y.subs{:}])
                    assign(reshape(OldArgument,[],1),reshape([Y(1).subs{:}],[],1));
                    varargout{1} = value(X);
                    return
                end
            end
            y = replace(X,OldArgument,[Y(1).subs{:}]);
            if isnumeric(y)
                varargout{1} = y;
                return
            end
            
        case '.'
            switch Y(1).subs
                case {'minimize','maximize'}
                    options = [];
                    constraints = [];                  
                    objective = varargin{1};
                    opsargs = {};
                    if length(Y)==2
                        if isequal(Y(2).type,'()')
                            for i = 1:length(Y(2).subs)
                                switch class(Y(2).subs{i})
                                    case {'lmi','constraint'}
                                        constraints = [constraints, Y(2).subs{i}];
                                    case 'struct'
                                        options = Y(2).subs{i};
                                    case {'double','char','gem','sgem'}
                                        opsargs{end+1} = Y(2).subs{i};
                                    otherwise
                                        error('Argument to minimize should be constraints or options');
                                end
                            end
                        else
                            error(['What do you mean with ' Y(2).type '?']);
                        end
                    end      
                    if length(opsargs)>0
                        if isempty(options)
                           options = sdpsettings(opsargs{:}); 
                        else
                            options = sdpsettings(options,opsargs{:});
                        end
                    end
                    if isequal(Y(1).subs,'minimize')
                        sol = solvesdp(constraints,objective,options);
                    else
                        sol = solvesdp(constraints,-objective,options);
                    end
                    varargout{1} = varargin{1};
                    varargout{2} = sol;
                    return
                case 'derivative'
                    try
                        m = model(varargin{1});
                        varargout{1} = m{1}.derivative;
                    catch
                        varargout{1} = 1;
                    end
                    return                
                otherwise
                    error(['Indexation  ''' Y.type Y.subs ''' not supported']) ;
            end
        otherwise
            error(['Indexation with ''' Y.type ''' not supported']) ;
    end
catch
    error(lasterr)
end
if isempty(y.lmi_variables)
    y = full(reshape(y.basis(:,1),y.dim(1),y.dim(2)));
else
    % Reset info about conic terms
    y.conicinfo = [0 0];
    y.extra.createTime = definecreationtime;
end
varargout{1} = y;

function X = subsref1d(X,ind1,Y)

% Get old and new size
n = X.dim(1);
m = X.dim(2);

% Convert to linear indecicies
if islogical(ind1)
    ind1 = double(find(ind1));
elseif ischar(ind1)
    X.dim(1) = n*m;
    X.dim(2) = 1;
    return;
elseif ~isnumeric(ind1)
    X = milpsubsref(X,Y);
    return
end

% Detect X(scalar)
if length(ind1) == 1 & ind1 <= n*m
    
    Z = X.basis(ind1,:);
    nnew = 1;
    mnew = 1;
    
else

    % What would the size be for a double
    dummy = reshape(X.basis(:,1),n,m);
    dummy = dummy(ind1);
    nnew = size(dummy,1);
    mnew = size(dummy,2);
    [nx,mx] = size(X.basis);
    
    if length(ind1) > 1
        try
            % row-based subsref for sparse objects can be very slow and
            % take a lot of memory, transposing a big sparse array as in
            %   Z = X.basis.';
            %   Z = Z(:,ind1);
            %   Z = Z.';
            % can also be quite slow for very large matrices, so we use
            % some custom code
            [ix,jx,sx] = find(X.basis);
            if length(ix) == length(unique(ix))
                % We never have two elements on the same row of X.basis
                [isNnz,loc] = ismember(ind1(:),ix);
                sel = loc(isNnz);
                ix = find(isNnz);
                jx = jx(sel);
                sx = sx(sel);
                Z = sparse(ix,jx,sx,numel(ind1),mx);
            else
                % We can have seveal elements on a given row of X.basis
                if numel(ind1) == length(unique(ind1))
                    % Every row of X.basis is selected at most once
                    [keep,loc] = ismember(ix,ind1(:));
                    ix = loc(keep);
                    jx = jx(keep);
                    sx = sx(keep);
                    Z = sparse(ix,jx,sx,numel(ind1),mx);
                else
                    % Complicated case, we use more generica code
                    Z = sparse(1:numel(ind1),ind1,1,numel(ind1),size(X,1))*X.basis;
                end
            end
            
            if false
                % For checking purpose
                Z0 = X.basis.';
                Z0 = Z0(:,ind1);
                Z0 = Z0.';
                assert(isequal(Z, Z0));
            end
        catch
            Z = X.basis(ind1,:);    
        end
    else
        Z = X.basis(ind1,:);
    end
end

% Find non-zero basematrices
nzZ = find(any(Z(:,2:end),1));
if ~isempty(nzZ)
    X.dim(1) = nnew;
    X.dim(2) = mnew;
    X.lmi_variables =  X.lmi_variables(nzZ);
    X.basis = Z(:,[1 1+nzZ]);
else
    bas = reshape(X.basis(:,1),n,m);
    X.dim(1) = nnew;
    X.dim(2) = mnew;
    X.lmi_variables = [];
    X.basis = reshape(bas(ind1),nnew*mnew,1);
end

function X = subsref2d(X,ind1,ind2,Y)

if isnumeric(ind1)
elseif ischar(ind1)
    ind1 = 1:X.dim(1);
elseif islogical(ind1)
    ind1 = double(find(ind1));
elseif ~isnumeric(ind1)
    X = milpsubsref(X,Y);
    return
end
if isnumeric(ind2)
elseif ischar(ind2)
    ind2 = 1:X.dim(2);
elseif islogical(ind2)
    ind2 = double(find(ind2));
elseif  ~isnumeric(ind2)    
    X = milpsubsref(X,Y);
    return
end

n = X.dim(1);
m = X.dim(2);
lind2 = length(ind2);
lind1 = length(ind1);
if lind2 == 1
    ind1_ext = ind1(:);
else
    ind1_ext = kron(ones(lind2,1),ind1(:));
end
if lind1 == 1
    ind2_ext = ind2(:);
else
    ind2_ext = kron(ind2(:),ones(lind1,1));
end

if any(ind1 > n) || any(ind2 > m)
    error('Index exceeds matrix dimensions.');
end
    
if lind1==1 && lind2==1
    if isequal(X.conicinfo,[-1 0])
        X.basis = [0 1];
        X.lmi_variables = X.lmi_variables(1)+ind1+(ind2-1)*n-1;
        X.dim = [1 1];
        X.conicinfo = [0 0];   
        return
    end
end

if prod(size(ind1_ext))==0 | prod(size(ind2_ext))==0
    linear_index = [];
else
    % Speed-up for some bizarre code with loads of indexing of vector
    if m==1 & ind2_ext==1
        linear_index = ind1_ext;
    elseif length(ind2_ext)==1 && length(ind1_ext)==1
        linear_index = ind1_ext + (ind2_ext-1)*n;
    else
        linear_index = sub2ind([n m],ind1_ext,ind2_ext);
    end
end
nnew = length(ind1);
mnew = length(ind2);

% Put all matrices in vectors and extract sub matrix
Z = X.basis(linear_index,:);
% Find non-zero basematrices
%nzZ = find(any(Z(:,2:end),1));
nzZ = find(any(Z,1))-1;
if numel(nzZ)>0
    if nzZ(1)==0
    nzZ = nzZ(2:end);
    end
end
if ~isempty(nzZ)
    X.dim(1) = nnew;
    X.dim(2) = mnew;
    X.lmi_variables =  X.lmi_variables(nzZ);
    X.basis = Z(:,[1 1+nzZ]);
else
    bas = reshape(X.basis(:,1),n,m);
    X.dim(1) = nnew;
    X.dim(2) = mnew;
    X.lmi_variables = [];
    X.basis = reshape(bas(linear_index),nnew*mnew,1);
end