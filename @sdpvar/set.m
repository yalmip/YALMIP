function F = set(varargin)
%SET OBSOLETE (but still used internally)

disp('SET has been considered obsolete for many years, and the time has come...');
disp('Update your code ');

error

if isa(varargin{1},'blkvar')
    varargin{1} = sdpvar(varargin{1});     
end

switch nargin
case 0
    F = lmi;
case 1
    F = lmi(varargin{1});
case 2
    F = lmi(varargin{1},varargin{2});
case 3
    F = lmi(varargin{1},varargin{1},varargin{3});
otherwise
end
    