function varargout = interp2_internal(varargin)
%INTERP2_INTERNAL (overloaded)

switch class(varargin{1})
    
    case 'double'
        if isequal(varargin{5},'graph') || isequal(varargin{5},'lp') || isequal(varargin{5},'milp') || isequal(varargin{5},'sos2')           
            varargout{1} = griddata(varargin{2},varargin{3},varargin{4},varargin{1}(1),varargin{1}(2));
            return
        end
        varargout{1} = interp2(varargin{2},varargin{3},varargin{4},varargin{1}(1),varargin{1}(2),varargin{5});
        
    case 'char'
        
        t = varargin{2};
        x = varargin{3}(1);
        y = varargin{3}(2);
        xi = varargin{4};
        yi = varargin{5};
        zi = varargin{6};
        method = varargin{7};
        
        % YALMIP requests a graph model
        if isequal(varargin{1},'graph')
            
            % Are we actually asking for a graph model?
            if strcmpi(method,'graph') || strcmpi(method,'lp')
                % Compute hyperplanes if possible
                [convexconcave,A,B] = isconvexmeshdata(xi,yi,zi);                
                if convexconcave == 1                
                    % Data convex, create epi-graph
                    F = [[x y t]*A + B >= 0, min(min(xi)) <= x <= max(max(xi)),min(min(yi)) <= y <= max(max(yi))];
                    mono = 'none';
                    def = 'none';
                    operator = struct('convexity','convex','monotonicity',mono,'definiteness',def,'model','graph');
                elseif convexconcave == -1       
                     % Data concave, create hypograph
                    F = [[x y t]*A + B <= 0, min(min(xi)) <= x <= max(max(xi)),min(min(yi)) <= y <= max(max(yi))];
                    mono = 'none';
                    def = 'none';                    
                    operator = struct('convexity','concave','monotonicity',mono,'definiteness',def,'model','graph');
                else
                    % Data is indefinite
                    F = [];
                    operator = [];
                end
            else
                F = [];
                operator = [];
            end
            varargout{1} = F;
            varargout{2} = operator;
            varargout{3} = varargin{3};
            
        else
            if isequal(varargin{end},'graph')
                % We want graph, YALMIP asks for exact representation
                varargout{1} = [];
                varargout{2} = [];
                varargout{3} = varargin{3};
            
            elseif strcmpi(method,'sos2') || strcmpi(method,'milp') || strcmpi(method,'lp')
                
                xi = varargin{4};
                yi = varargin{5};
                zi = varargin{6};
                [N,M] = size(xi);
                xy = varargin{3};
                x = xy(1);
                y = xy(2);
                z = varargin{2};
                
                Lambda = sdpvar(N,M,'full');
                xinterp = sum(sum(xi.*Lambda));
                yinterp = sum(sum(yi.*Lambda));
                zinterp = sum(sum(zi.*Lambda));
                colsos = sdpvar(M,1,'full');
                rowsos = sdpvar(N,1,'full');
                F = [sum(sum(Lambda))==1, Lambda >= 0];
                F = [F, colsos == sum(Lambda,1)'];
                F = [F, rowsos == sum(Lambda,2)];
                F = [F, 0 <= colsos <= 1, 0 <= rowsos <= 1];
                F = [F, sos2(colsos), sos2(rowsos)];
                F = [F, z == sum(sum(Lambda.*zi)),
                    x == sum(sum(Lambda.*xi)),
                    y == sum(sum(Lambda.*yi))];
                
                Elements = reshape(1:N*M,N,M);
                for i = -(N+M):N+M
                    j = diag(fliplr(Elements),i);
                    if length(j) >= 3
                        F = [F, sos2(Lambda(j))];
                    end
                end
                operator = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
                varargout{1} = F;
                varargout{2} = operator;
                varargout{3} = varargin{3};
                
            else
                L = [min(varargin{4}(:));min(varargin{5}(:))];
                U = [max(varargin{4}(:));max(varargin{5}(:))];
                operator = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');
                operator.bounds = @bounds;
                varargout{1} = [L <= varargin{3} <= U];
                varargout{2} = operator;
                varargout{3} = varargin{3};
            end
        end
    otherwise
        error('SDPVAR/INTERP2_INTERNAL called with strange argument');
end

function [L,U] = bounds(xL,xU,varargin)

xv = varargin{1};
yv = varargin{2};
zv = varargin{3};
i = 1;
lowindex = zeros(2,1);
highindex = zeros(2,1);
if xL(1) <= xv(1)
    lowindex(i) = 1;
else
    temp = find(xL(1) > xv(1,:));lowindex(i) = max(1,max(temp));
end
if xU(1) >= xv(end)
    highindex(i) = size(xv,2);
else
    temp = find(xv(1,:) < xU(1));highindex(i) = min(size(xv,2),max(temp)+1);
end
i = 2;
if xL(2) <= yv(1)
    lowindex(i) = 1;
else
    temp = find(xL(2) > yv(:,1));lowindex(i) = max(1,max(temp));
end
if xU(2) >= yv(end)
    highindex(i) = size(yv,1);
else
    temp = find(yv(:,1) < xU(2));highindex(i) = min(size(yv,1),max(temp)+1);
end

if isequal(varargin{4},'linear')
    Z = zv(lowindex(1):highindex(1),lowindex(2):highindex(2));
    L = min(Z(:));
    U = max(Z(:));
else
    dx = (xv(1,highindex(1))-xv(1,lowindex(1)))/(size(xv,2)*25);
    dy = (yv(highindex(2),1)-yv(lowindex(2),1))/(size(yv,1)*25);
    [X,Y] = meshgrid(xv(1,lowindex(1)):dx:xv(1,highindex(1)),yv(lowindex(2),1):dy:yv(highindex(2),1));
    Z = interp2(xv,yv,zv,X,Y,varargin{4});
    % To account for finite grid, we add a precision dependent margin
    L = min(Z(:))-dx-dy;
    U = max(Z(:))+dx+dy;
end