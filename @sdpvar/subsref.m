function varargout = subsref(varargin)
%SUBSREF (overloaded)

% Author Johan Löfberg
% $Id: subsref.m,v 1.23 2009-09-28 07:36:25 joloef Exp $


% Stupid first slice call (supported by MATLAB)
% x = sdpvar(2);x(1,:,:)
if length(varargin{2})==1
    if  length(varargin{2}.subs) > 2 && isequal(varargin{2}.type,'()')
        i = 3;
        ok = 1;
        while ok && (i <= length(varargin{2}.subs))
            ok = ok && (isequal(varargin{2}.subs{i},1) || isequal(varargin{2}.subs{i},':'));
            i = i + 1;
        end
        if ok
            varargin{2}.subs = {varargin{2}.subs{1:2}};
        else
            error('??? Index exceeds matrix dimensions.');
        end
    end
end

if (isequal(varargin{2}.type,'()') && ((isa(varargin{2}.subs{1},'sdpvar')) || (length(varargin{2}.subs)==2 && isa(varargin{2}.subs{2},'sdpvar'))))
    % *****************************************
    % Experimental code for varaiable indicies
    % *****************************************
    varargout{1} = milpsubsref(varargin{:});
    return
else
    X = flush(varargin{1});
    Y = varargin{2};
end

try
    switch Y(1).type
        case '()'
            if  isa(Y(1).subs{1},'constraint')
                error('Conditional indexing not supported.');
            end
            % Check for simple cases to speed things up (yes, ugly but we all want speed don't we!)
            switch size(Y(1).subs,2)
                case 1
                    if isa(Y(1).subs{1},'sdpvar')
                        varargout{1} = yalmip('addextendedvariable',mfilename,varargin{:});
                        return
                    else
                        y = subsref1d(X,Y(1).subs{1});
                    end
                case 2
                    y = subsref2d(X,Y.subs{1},Y(1).subs{2});
                otherwise
                    error('Indexation error.');
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
                    XX = double(X);
                    varargout{1} = double(X);
                    return
                end

%                     vars = getvariables(X(i));
%                     nonlinearModel = yalmip('extstruct',vars);
%                     for j = 1:length(nonlinearModel)
%                         if ~((isequal(nonlinearModel{j}.fcn,'pwa_yalmip') | isequal(nonlinearModel{j}.fcn,'pwq_yalmip'))& isa(Y.subs{:},'double'))
%                             mpt_solution = 0;
%                         end
%                     end
%                     if mpt_solution
%                         assign(nonlinearModel{1}.arg{2},Y.subs{:});
%                         XX = double(X);
%                         varargout{1} = [varargout{1};XX(i)];
%                     end
%                 end
%                 if mpt_solution
%                     return;
%                 end
            end
            vars = getvariables(X);
            if (length(vars) == 1) & ismembc(vars,yalmip('extvariables'))
                nonlinearModel = yalmip('extstruct',vars);
%                 if (isequal(nonlinearModel.fcn,'pwa_yalmip') | isequal(nonlinearModel.fcn,'pwq_yalmip'))& isa(Y.subs{:},'double')
%                     assign(nonlinearModel.arg{2},Y.subs{:});
%                     varargout{1} = double(X);
%                     return
%                 end
                OldArgument = [];
                for i = 1:length(nonlinearModel.arg)
                    if isa(nonlinearModel.arg{i},'sdpvar')
                        OldArgument = [OldArgument;  nonlinearModel.arg{i}];
                    end
                end
                if isa([Y.subs{:}],'double')
                    assign(reshape(OldArgument,[],1),reshape([Y(1).subs{:}],[],1));
                    varargout{1} = double(X);
                    return
                end
            end
            y = replace(X,OldArgument,[Y(1).subs{:}]);
            if isa(y,'double')
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
                                    case {'double','char'}
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
end
varargout{1} = y;

function X = subsref1d(X,ind1)

% Get old and new size
n = X.dim(1);
m = X.dim(2);

% Convert to linear indecicies
if islogical(ind1)
    ind1 = double(find(ind1));
end

% Ugly hack handle detect X(:)
%pickall = 0;
if ischar(ind1)
    X.dim(1) = n*m;
    X.dim(2) = 1;
    return;
end

% Detect X(scalar)
if length(ind1) == 1 & ind1 <= n*m
    
    Z = X.basis.';
    Z = Z(:,ind1);
    Z = Z.';
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
            Z = X.basis.';
            Z = Z(:,ind1);
            Z = Z.';
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

function X = subsref2d(X,ind1,ind2)

if ischar(ind1)
    ind1 = 1:X.dim(1);
end
if ischar(ind2)
    ind2 = 1:X.dim(2);
end

% Convert to linear indecicies
if islogical(ind1)
    ind1 = double(find(ind1));
end

% Convert to linear indecicies
if islogical(ind2)
    ind2 = double(find(ind2));
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

if prod(size(ind1_ext))==0 | prod(size(ind2_ext))==0
    linear_index = [];
else
    % Speed-up for some bizarre code with loads of indexing of vector
    if m==1 & ind2_ext==1
        linear_index = ind1_ext;
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