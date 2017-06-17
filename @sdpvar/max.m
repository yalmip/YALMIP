function [y,loc] = max(varargin)
%MAX (overloaded)
%
% t = max(X)
% t = max(X,Y)
% t = max(X,[],DIM)
%
% Creates an internal structure relating the variable t with convex
% operator max(X).
%
% The variable t is primarily meant to be used in convexity preserving
% operations such as t<=..., minimize t etc.
%
% If the variable is used in a non-convexity preserving operation, such as
% t>=0, a mixed integer model will be derived.
%
% See built-in MAX for syntax.

% To simplify code flow, code for different #output/inputs
if nargout == 2
    [y,loc] = max_with_loc(varargin{:});
    return
end
        
switch nargin
    case 1
        % Four cases:
        % 1. One scalar input, return same as output
        % 2. A vector input should give scalar output
        % 3. Matrix input returns vector output
        % 4. User wants location index
                      
        X = varargin{1};
        
        if max(size(X))==1
            y = X;
            return
        elseif min(size(X))==1
            X = removeInf(X);
            if isnumeric(X)
                y = max(X);
            elseif length(X) == 1
                y = X;
            else
                
                y = yalmip('define','max_internal',X);
                % Some special code to ensure max(x) when x is a simple
                % binary vector yields a binary graph variable. This will
                % simplify some models
                reDeclareForBinaryMax(y,X);
            end
            return
        else
            % This is just short-hand for general command
            y = max(X,[],1);
        end
        
    case 2
        
        X = varargin{1};
        Y = varargin{2};
        [nx,mx] = size(X);
        [ny,my] = size(Y);
        if ~((nx*mx==1) | (ny*my==1))
            % No scalar, so they have to match
            if ~((nx==ny) & (mx==my))
                error('Array dimensions must match.');
            end
        end
        
        % Convert to compatible matrices
        if nx*mx==1
            X = X*ones(ny,my);
            nx = ny;
            mx = my;
        elseif ny*my == 1
            Y = Y*ones(nx,mx);
            ny = nx;
            my = mx;
        end
        
        % Ok, done with error checks etc.
        Z = [reshape(X,1,[]);reshape(Y,1,[])];
        y = yalmip('define','max_internal',Z);
        reDeclareForBinaryMax(y,Z);
        y = reshape(y,nx,mx);
        
    case 3
        
        X = varargin{1};
        Y = varargin{2};
        DIM = varargin{3};
        
        if ~(isa(X,'sdpvar') & isempty(Y))
            error('MAX with two matrices to compare and a working dimension is not supported.');
        end
        
        if ~isa(DIM,'double')
            error('Dimension argument must be 1 or 2.');
        end
        
        if ~(length(DIM)==1)
            error('Dimension argument must be 1 or 2.');
        end
        
        if ~(DIM==1 | DIM==2)
            error('Dimension argument must be 1 or 2.');
        end
        
        if DIM==1
            % Create one extended variable per column
            y = [];
            for i = 1:size(X,2)
                inparg = extsubsref(X,1:size(X,1),i);
                if isa(inparg,'sdpvar')
                    inparg = removeInf(inparg);
                    if  isnumeric(inparg)
                        y = [y max(inparg)];
                    elseif length(inparg) == 1
                        y = [y max(inparg)];
                    else
                        z = yalmip('define','max_internal',inparg);
                        y = [y z];
                        % Some special code to ensure max(x) when x is a simple
                        % binary vector yields a binary graph variable. This will
                        % improve some models
                        reDeclareForBinaryMax(z,inparg);
                    end
                else
                    y = [y max(inparg)];
                end
            end
        else
            % Re-use code recursively
            y = max(X',[],1)';
        end
        
    otherwise
        error('Too many input arguments.');
end

function X = removeInf(X)
Xbase = getbase(X);
infs = find( isinf(Xbase(:,1)) & (Xbase(:,1)<0));
if ~isempty(infs)
    X.basis(infs,:) = [];
    if X.dim(1)>X.dim(2)
        X.dim(1) = X.dim(1) - length(infs);
    else
        X.dim(2) = X.dim(2) - length(infs);
    end
end
infs = find(isinf(Xbase(:,1)) & (Xbase(:,1)>0));
if ~isempty(infs)
    X = inf;
end