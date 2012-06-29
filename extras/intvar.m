function sys = intvar(varargin)
%INTVAR Create symbolic integer variable 
%   
%   INTVAR works exactly as SDPVAR, with the only difference that
%   the elements in the variable automatically will be constrained 
%   to be integer when an optimization problem is solved.
%
%   See also SDPVAR, BINVAR, INTEGER, BINARY

% Author Johan Löfberg
% $Id: intvar.m,v 1.7 2006-05-17 13:22:29 joloef Exp $

all_variable_names = 1;
i = 1;
while all_variable_names & i<=nargin
    all_variable_names = all_variable_names & isvarname(varargin{i});
    i = i + 1;
end
if all_variable_names
    for k = 1:nargin
        varname = varargin{k};
        assignin('caller',varname,intvar(1,1));
    end
else
    sys = sdpvar(varargin{:}); 
    if isa(sys,'cell')
        for i = 1:length(sys)
            yalmip('setintvariables',[yalmip('intvariables') getvariables(sys{i})]);
        end
    else
        yalmip('setintvariables',[yalmip('intvariables') getvariables(sys)]);
    end
end