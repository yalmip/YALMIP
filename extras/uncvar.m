function sys = uncvar(varargin)
%UNCVAR Create symbolic uncertain variable 
%   
%   UNCVAR works exactly as SDPVAR, [FIX...]
%
%   See also SDPVAR

all_variable_names = 1;
i = 1;
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
    switch nargin
        case 1
            sys = sdpvar(varargin{1});
        case 2
            sys = sdpvar(varargin{1},varargin{2});
        case 3
            sys = sdpvar(varargin{1},varargin{2},varargin{3});
        case 4
            sys = sdpvar(varargin{1},varargin{2},varargin{3},varargin{4});
        otherwise
            error('Wrong number of input arguments. See help-text for sdpvar')
    end
    yalmip('setuncvariables',[yalmip('uncvariables') getvariables(sys)]);
end