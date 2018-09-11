function varargout = interp1(varargin)
%INTERP1 (overloaded)

switch class(varargin{3})

    case 'sdpvar'   
                
        if nargin < 4
            varargin{4} = 'linear';
        end
        
        if nargin < 5
            varargin{5} = [];
        end        
        
        if isa(varargin{2},'function_handle')
            varargin{2} = varargin{2}(varargin{1});
        end
        
        if numel(varargin{3}) > 1
            d = size(varargin{3});
            
            if min(d) > 1
                % Matrix case, vectorize
                varargout{1} = interp1(varargin{1},varargin{2},reshape(varargin{3},[],1),varargin{4},varargin{5});
                varargout{1} = reshape(varargout{1},d);
                return
            end
            
            if any(any(isnan(varargin{1})))
                error('Interpolation grid contains NaNs');
            end
            if any(any(isnan(varargin{2})))
                error('Interpolation data contains NaNs');
            end
                     
            % If only one data sequence, place in row vector
            if min(size(varargin{1})) == 1
                varargin{1} = repmat(varargin{1}(:)',max(d),1);
            end
            if min(size(varargin{2})) == 1
                varargin{2} = repmat(varargin{2}(:)',max(d),1);
            end           
            
            out = [];
            varargin{3} = reshape(varargin{3},1,max(d));
            for i = 1:max(d)
                    Y.type = '()';
                    Y.subs{1} = 1;
                    Y.subs{2} = i;
                    xij = subsref(varargin{3},Y);                    
                    out = [out;interp1(varargin{1}(i,:),varargin{2}(i,:),xij,varargin{4},varargin{5})];                                
            end
            if d(1)==1
                out = out';
            end
            varargout{1} = out;
            return
        else
            if min(size(varargin{1})) > 1 || min(size(varargin{2})) > 1
                error('Data should be vectors for scalar interpolant');
            end            
            % Column vector assumed format
            varargin{1} = reshape(varargin{1},[],1);
            varargin{2} = reshape(varargin{2},[],1);            
        end
        
        if ~isa(varargin{1},'double') || ~isa(varargin{2},'double')
            error('First 2 arguments in interp1 should be approximation data');
        end
        
        %if any(diff(varargin{1})<0)
        %    error('First arguments in interp1 should be monotonically increasing');
        %end
        
        if isequal(varargin{4},'graph')
            if  ~isconvexdata(varargin{1},varargin{2}) && ~isconvexdata(varargin{1},-varargin{2})
                error('Data has to be convex or concave for graph approximant');
            end
        end
        
        % Reorder arguments to make sure sdpvar is first argument.
        % Use local version of interp to deal with this
        % Arguments xi yi x method -> x xi yi method
        varargout{1} = yalmip('define','interp1_internal',varargin{[3 1 2 4 5]});   
     
    case 'double'
        if isa(varargin{1},'double') && isa(varargin{2},'sdpvar')             
            % User is trying to do nonlinear programming where the knot
            % function values in an interpolant are free variables
            zi = varargin{3};
            [n,m] = size(zi);
            if n*m == 1
                varargout{1} = yalmip('define','interp1_nonlinear',varargin{1},varargin{2},varargin{3},varargin{4});
            else
                out = [];
                zi = reshape(zi,1,n*m);
                for i = 1:n*m
                    Y.type = '()';
                    Y.subs{1} = 1;
                    Y.subs{2} = i;
                    zij = subsref(zi,Y);
                    out = [out  yalmip('define','interp1_nonlinear',varargin{1},varargin{2},zij,varargin{4})];
                end
                varargout{1} = reshape(out,n,m);
            end

                            
        elseif isa(varargin{1},'sdpvar') && isa(varargin{2},'sdpvar') 
            % User is trying to do nonlinear programming where the knot
            % function values and locations in an interpolant are free variables
            
            zi = varargin{3};
            [n,m] = size(zi);
            if n*m == 1
                varargout{1} = yalmip('define','interp1_nonlinear',[varargin{1},varargin{2}],[],varargin{3},varargin{4});
            else
                out = [];
                zi = reshape(zi,1,n*m);
                for i = 1:n*m
                    Y.type = '()';
                    Y.subs{1} = 1;
                    Y.subs{2} = i;
                    zij = subsref(zi,Y);
                    out = [out   yalmip('define','interp1_nonlinear',[varargin{1},varargin{2}],[],zij,varargin{4})];
                end
                varargout{1} = reshape(out,n,m);
            end            
        else
            error('SDPVAR/INTERP1 called with strange arguments!');    
        end
        
    otherwise
        error('SDPVAR/INTERP1 called with strange argument!');
end

function isconvex = isconvexdata(xi,yi)
finterp = (yi(1:end-2) + yi(3:end))/2;
if all(finterp >= yi(2:end-1))
    isconvex = 1;
else
    isconvex = 0;
end