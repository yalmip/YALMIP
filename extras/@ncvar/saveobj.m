function out = saveobj(obj)
%SAVEOBJ (overloaded)

% Author Johan Löfberg 
% $Id: saveobj.m,v 1.1 2006-08-10 18:00:22 joloef Exp $   

% We have to save the persistent variables in the SDPVAR class
obj.savedata.internalsdpvarstate = yalmip('getinternalsdpvarstate');
obj.savedata.version = yalmip('version');
out = obj;
