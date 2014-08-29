function varargout = subsref(varargin)
%SUBSREF (overloaded)

% Stupid first slice call (supported by MATLAB)
if  length(varargin{2}.subs) > 2
    i = 3;
    ok = 1;
    while ok & (i <= length(varargin{2}.subs))
        ok = ok & (isequal(varargin{2}.subs{i},1) | isequal(varargin{2}.subs{i},':'));
        i = i + 1;
    end
    if ok
        varargin{2}.subs = {varargin{2}.subs{1:2}};
    else
        error('??? Index exceeds matrix dimensions.');
    end
    
end

if (isa(varargin{2}.subs{1},'sdpvar')) | (length(varargin{2}.subs)==2 & isa(varargin{2}.subs{2},'sdpvar'))
    % *****************************************
    % Experimental code for varaiable indicies   
    % *****************************************
    varargout{1} = milpsubsref(varargin{:});
    return
else
    X = varargin{1};
    Y = varargin{2};
end

try
    switch Y.type
        case '()'
            % Check for simple cases to speed things up (yes, ugly but we all want speed don't we!)
            switch size(Y.subs,2)
                case 1
                    if isa(Y.subs{1},'sdpvar')
                        varargout{1} = yalmip('addextendedvariable',mfilename,varargin{:});
                        return
                    else
                        y = subsref1d(X,Y.subs{1});
                    end
                case 2
                    y = subsref2d(X,Y.subs{1},Y.subs{2});
                otherwise
                    error('Indexation error.');
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

% What would the size be for a double
dummy = reshape(X.basis(:,1),n,m);
dummy = dummy(ind1);
nnew = size(dummy,1);
mnew = size(dummy,2);

[nx,mx] = size(X.basis);
% Sparse row-based subsref can be EEEEEEEXTREMELY SLOW IN SOME CASES
% FIX : Smarter approach?
% 
% try
%     [ix,jx,sx] = find(X.basis);
%     [keep ,loc] = ismember(ix,ind1);keep = find(keep);
%     ix = loc(keep);%ix = loc(ix);
%     jx = jx(keep);
%     sx = sx(keep);
%     Z = sparse(ix,jx,sx,length(ind1),mx);
% catch
%     9
% end

if length(ind1) > 1
    Z = X.basis.';
    Z = Z(:,ind1);
    Z = Z.';
else
    Z = X.basis(ind1,:);
end

% if ~isequal(Z,Z2)
%     error
% end

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
%nzZ = find(any(Z,1));
% % A bit messy code to speed up things
% if isempty(nzZ)
%     bas = reshape(X.basis(:,1),n,m);
% 	X.dim(1) = nnew;
% 	X.dim(2) = mnew;
% 	X.lmi_variables = [];
% 	X.basis = reshape(bas(ind1),nnew*mnew,1);    
% else
%     if nzZ(1) == 1
%         
%     else
%     end
% end

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
    ind1_ext = kron(repmat(1,lind2,1),ind1(:));
end
if lind1 == 1
    ind2_ext = ind2(:);
else
    ind2_ext = kron(ind2(:),repmat(1,lind1,1));
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
    X.basis = reshape(bas(linear_index),nnew*mnew,1);
end