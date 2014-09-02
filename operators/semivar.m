function sys = semivar(varargin)
%SEMIVAR Create semicontinuous variable
%
%   SEMIVAR works exactly as SDPVAR, with the only difference that
%   the elements in the variable automatically will be constrained
%   to be semicontinous (x = 0 or [L <= x <= U])
%
%   Thwe lower and upper bounds are defined separately.
%
%   See also SDPVAR, BINVAR, INTVAR, BINARY, INTEGER

all_variable_names = 1;
i = 1;
while all_variable_names & i<=nargin
    all_variable_names = all_variable_names & isvarname(varargin{i});
    i = i + 1;
end
if all_variable_names
    for k = 1:nargin
        varname = varargin{k};
        assignin('caller',varname,semivar(1,1));
    end
else
    sys = sdpvar(varargin{:}); 
    if isa(sys,'cell')
        for i = 1:length(sys)
            yalmip('setsemicontvariables',[yalmip('semicontvariables') getvariables(sys{i})]);
        end
    else
        yalmip('setsemicontvariables',[yalmip('semicontvariables') getvariables(sys)]);
    end
end