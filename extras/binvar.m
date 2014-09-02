function sys = binvar(varargin)
%BINVAR Create symbolic binary variable
%
%   BINVAR works exactly as SDPVAR, with the only difference that
%   the elements in the variable automatically will be constrained
%   to be binary (0/1) when an optimization problem is solved.
%
%   See also INTVAR, SDPVAR, BINARY, INTEGER

all_variable_names = 1;
i = 1;

% check for command line syntax
while all_variable_names & i<=nargin
    all_variable_names = all_variable_names & isvarname(varargin{i});
    i = i + 1;
end

if all_variable_names
    for k = 1:nargin
        varname = varargin{k};
        assignin('caller',varname,binvar(1,1));
    end
else
    sys = sdpvar(varargin{:});
    if isa(sys,'cell')
        for i = 1:length(sys)
            yalmip('setbinvariables',[yalmip('binvariables') getvariables(sys{i})]);
        end
    else
        yalmip('setbinvariables',[yalmip('binvariables') getvariables(sys)]);
    end
end