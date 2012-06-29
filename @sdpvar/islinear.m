function sys=islinear(x)
%ISLINEAR Check if variable is linear       
%
% p = islinear(X)
%
% X : SDPVAR object
% p : boolean 0/1

% Author Johan Löfberg 
% $Id: islinear.m,v 1.2 2004-07-02 08:17:29 johanl Exp $  

try
    sys = is(x,'linear');
catch
    disp('I have messed up something internally in YALMIP. Sorry...')
    disp('Type yalmip(''clear'') and re-define the problem.');
    disp(' ')
    disp('This problem typically occurs if you have SDPVAR variables in');
    disp('work-space, and call a function where you are using ')
    disp('the command yalmip(''clear''). When control is returned to');
    disp('the work-space, things are messed up...');
    error('Internal problems in YALMIP');
end