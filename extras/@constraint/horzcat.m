function F = horzcat(varargin)

F = [];
for i=1:1:nargin
    if isa(varargin{i},'double') & ~isempty(varargin{i})
        warning('One of the constraints evaluates to a DOUBLE variable');
    elseif isa(varargin{i},'logical')
        if all(varargin{i}==1)
            % warning('One of the constraints evaluates to a LOGICAL variable');
        else
            error('One of the constraints evaluates to a FALSE LOGICAL variable. Your model is infeasible');
        end
    elseif isa(varargin{i},'optproblem')
        F = [varargin{i},F,varargin{i+1:end}];
        return;
    else
        H = lmi(varargin{i});
        F = F+H;
    end
end