function varargout = interp1_internal(varargin)
%INTERP1_INTERNAL (overloaded)

switch class(varargin{1})
    
    case 'double'
        if isequal(varargin{4},'graph') || isequal(varargin{4},'lp') || isequal(varargin{4},'milp') || isequal(varargin{4},'sos2')
            varargout{1} = interp1(varargin{2},varargin{3},varargin{1},'linear');
        else
            varargout{1} = interp1(varargin{2},varargin{3},varargin{1},varargin{4});
        end
        
    case 'char'
        
        t = varargin{2};
        x = varargin{3};
        xi = varargin{4};
        yi = varargin{5};
        method = varargin{6};
        dyi = varargin{7};
        
        if strcmpi(varargin{1},'graph')
            % Convexity propagation wants epi-graph models, so is that what
            % we have used for modelling?
            if strcmpi(method,'envelope')
                   [xi,yi] = convhullprune(xi,yi);
                   [Model,mono,def] = convexGraph(xi,yi,x,t,dyi);  
                   operator = struct('convexity','convex','monotonicity',mono,'definiteness',def,'model','graph');                    
            elseif strcmpi(method,'graph') || strcmpi(method,'lp')
                if isconvexdata(xi,yi)
                    % Convex case, create epi-graph
                    [Model,mono,def] = convexGraph(xi,yi,x,t,dyi);                                          
                    operator = struct('convexity','convex','monotonicity',mono,'definiteness',def,'model','graph');
                elseif isconvexdata(xi,-yi)
                    [Model,mono,def] = concaveGraph(xi,yi,x,t,dyi);                      
                    operator = struct('convexity','concave','monotonicity',mono,'definiteness',def,'model','graph');
                else
                    Model = [];
                    operator = [];
                end
            else
                Model = [];
                operator = [];
            end
        else
            if strcmpi(method,'sos2') || strcmpi(method,'milp') || strcmpi(method,'lp')
                lambda = sdpvar(length(xi),1);
                Model = [sos2(lambda), x == lambda'*xi(:), t == lambda'*yi(:),lambda>=0, sum(lambda)==1];
                [mono,def] = classifyData(yi);
                operator = struct('convexity','none','monotonicity',mono,'definiteness',def,'model','integer');
            elseif strcmpi(method,'graph')
                Model = [];
                operator = [];
            else
                Model = [min(xi) <= varargin{3} <= max(xi)];
                operator = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');
                operator.bounds = @bounds;
            end
        end
        
        varargout{1} = Model;
        varargout{2} = operator;
        varargout{3} = varargin{3};
        
    otherwise
        error('SDPVAR/INTERP1_INTERNAL called with strange argument');
end

function [L,U] = bounds(xL,xU,varargin)

xv = varargin{1};
yv = varargin{2};
if xL <= xv(1)
    index1 = 1;
else
    index1 = find(xL > xv);index1 = max(1,max(index1));
end
if xU >= xv(end)
    index2 = length(xv);
else
    index2 = find(xv < xU);index2 = min(length(xv),max(index2)+1);
end

if isequal(varargin{3},'linear')
    Z = yv(index1:index2);
    L = min(yv(index1:index2));
    U = max(yv(index1:index2));
else
    N = ceil((xv(index2)-xv(index1))/(mean(diff(xv))/100));
    z = linspace(xv(index1),xv(index2),N);
    yz = interp1(xv,yv,z,varargin{3});
    % To account for finite grid, we add a precision dependent margin
    dz = (z(2)-z(1));
    L = min(yz)-dz;
    U = max(yz)+dz;
end

function isconvex = isconvexdata(xi,yi)
finterp = (yi(1:end-2) + yi(3:end))/2;
if all(finterp >= yi(2:end-1))
    isconvex = 1;
else
    isconvex = 0;
end

function [Model,mono,def] = convexGraph(xi,yi,x,t,dyi)


Model = [min(xi) <= x <= max(xi)];
if isempty(dyi)
    grads = diff(yi)./diff(xi);
    Model = [Model, yi(1:end-1) + grads.*(x - xi(1:end-1)) <= t];
else
    grads = dyi;
    Model = [Model, yi(:) + grads(:).*(x - xi(:)) <= t];
end
[mono,def] = classifyData(yi);

function [Model,mono,def] = concaveGraph(xi,yi,x,t,dyi)

Model = [min(xi) <= x <= max(xi)];
if isempty(dyi)
    grads = diff(yi)./diff(xi);
    Model = [Model, yi(1:end-1) + grads.*(x - xi(1:end-1)) >= t];
else
    grads = dyi;
    Model = [Model, yi(:) + grads(:).*(x - xi(:)) >= t];
end

[mono,def] = classifyData(yi);

function [xout,yout] = convhullprune(xin,yin)
% Prunes data down to the lower convex envelope

% First prune by convex hull
cj = convhulln([xin(:) yin(:)]);
cj = sort(unique(cj));
x = xin(cj);
y = yin(cj);

% This can probably be done in linear complexity, I am just lazy here
keep = ones(length(x),1);
for i = 1:length(x)-1
    if keep(i)
        for j = length(x):-1:3
            if keep(j)
                grad = (y(j)-y(i))/(x(j)-x(i));
                k = i+1:j-1;
                fail =  y(k) > y(i) + grad*(x(k)-x(i));
                keep(k(fail)) = 0;
            end
        end
    end
end
keep = cj(find(keep));
xout = xin(keep);
yout = yin(keep);


function [mono,def] = classifyData(yi);
if all(diff(yi)) >= 0
    mono = 'increasing';
elseif all(diff(yi) <= 0)
    mono = 'decreasing';
else
    mono = 'none';
end
if all(yi>=0)
    def = 'positive';
elseif all(yi <= 0)
    def = 'negative';
else
    def = 0;
end
