function deg=degree(p,y,flag,vector)
%DEGREE Polynomial degree
%
% DEG = DEGREE(p,x,flag,vector)
%
% p      : SDPVAR object.
% x      : Degree w.r.t linear SDPVAR objects.
% flag   : 'max', 'min'. Default 'max'
% vector : If vector = 1, returns degree of each element in p
%
% Examples
% x1 = sdpvar(1,1);x2 = sdpvar(1,1);
% p = [x1;x1*x2+x2^2];
%
% degree(p) returns 2
%
% degree(p,x1) returns 1
%
% degree(p,[x1 x2]) returns [1 2]
%
% degree(p,[x1 x2],[],1) returns [1 0;1 2]
%
% degree(p,[],1) returns [1;2]

if nargin == 3
    if isa(flag,'double')       
        % Old syntax
        deg = degree(p,y,'max',flag);
        return;
    end
   
end

if isnumeric(p)
    if nargin==1
        deg = 0;
    else
        deg = zeros(1,length(y));
    end
    return
end

if nargin<2 || isempty(y)
    y = recover(depends(p));   
end

if nargin < 3
    flag = 'max';
end

if nargin < 4
    vector = 0;
end

if vector == 0
    exponent_p = exponents(p,y);
    switch nargin
        case 1            
            degrees = sum(exponent_p,2);
            deg = full(max(degrees));            
        case {2,3}
            switch flag
                case 'max'
                    deg = full(max(exponent_p,[],1));
                case 'min'                    
                    deg = full(min(exponent_p,[],1));
                otherwise
            end
        otherwise
            error('Too many arguments. Wadda ya mean?')
    end
else
    p = p(:);
    if isempty(y)
        yy = recover(depends(p));
    else
        yy = y;
    end
    
    deg = zeros(length(p),length(depends(yy)));
    for i = 1:length(p)
        z.type = '()';
        z.subs{1} = i;
        pi = subsref(p,z);
        if isnumeric(pi)
            deg(i,:) = 0;
        else
            exponent_p = exponents(pi,yy);
            switch nargin
                case 1
                    deg(i,:) = full(max(sum(exponent_p,2)));
                case {2,3,4}
                    switch flag
                        case 'max'
                            deg(i,:) = full(max(exponent_p,[],1));
                        case 'min'                           
                            deg(i,:) = full(min(exponent_p,[],1));
                        otherwise
                            error('Do not understand the flag')
                    end
                otherwise
                    error('Too many arguments. Wadda ya mean?')
            end
        end
    end
    if isempty(y)
        deg = sum(deg,2);
    end
end