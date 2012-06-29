function F = vertcat(varargin)

F = set(varargin{1});
for i=2:1:nargin
    if isa(varargin{i},'double') 
       warning('One of the constraints evaluates to a DOUBLE variable');
    elseif isa(varargin{i},'logical') 
         warning('One of the constraints evaluates to a LOGICAL variable');
    else
        F=F+set(varargin{i});
    end
end