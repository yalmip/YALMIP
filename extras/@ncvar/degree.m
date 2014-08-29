function deg=degree(p,y,e)
%DEGREE Polynomial degree
%
% DEG = DEGREE(p,x,e)
%
% p : SDPVAR object.
% x : Degree w.r.t linear SDPVAR objects, can be [].


% e : If e=1, returns degree of each element in p
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
% degree(p,[x1 x2],1) returns [1 0;1 2]
%
% degree(p,[],1) returns [1;3]

if isa(p,'double')
    if nargin==1
        deg = 0;
    else
        deg = zeros(1,length(y));
    end
    return
end

if nargin<2
    y = recover(depends(p));
end

if nargin<3 | (nargin==3 & e==0)
    exponent_p = exponents(p,y);
    switch nargin
        case 1
            deg = full(max(sum(exponent_p,2)));
        case {2,3}
            deg = full(max(exponent_p,[],1));
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
    
    for i = 1:length(p)
        z.type = '()';
        z.subs{1} = i;
        exponent_p = exponents(subsref(p,z),yy);       
        switch nargin
            case 1
                deg(i,:) = full(max(sum(exponent_p,2)));
            case {2,3}
                deg(i,:) = full(max(exponent_p,[],1));
            otherwise
               error('Too many arguments. Wadda ya mean?')
        end
    end
    if isempty(y)
        deg = sum(deg,2);
    end
end