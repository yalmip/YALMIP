function F = vertcat(varargin)

F = set(varargin{1});
for i=2:1:nargin
    if isa(varargin{i},'double')
        warning('One of the constraints evaluates to a DOUBLE variable');
    elseif isa(varargin{i},'logical')
        if all(varargin{i}==1)
            warning('One of the constraints evaluates to a LOGICAL variable');
        else
            warning('One of the constraints evaluates to a FALSE LOGICAL variable');
        end
    elseif isa(varargin{i},'optproblem')
        F = [varargin{i},F,varargin{i+1:end}];
        return;
    else
        F=F+set(varargin{i});
    end
end