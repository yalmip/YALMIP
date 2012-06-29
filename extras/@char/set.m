function F = set(varargin)
switch nargin
case 0
    F = lmi;
case 1
    % Yihaa! 
    F = evalin('caller',['lmi(' '''' strrep(varargin{1},'''','''''') '''' ')']);
case 2
    F = evalin('caller',['lmi(' strrep(varargin{1},'''','''''') ','  '''' varargin{2} '''' ')']);
otherwise
end
    