function out = saveobj(obj)
%SAVEOBJ Save filter for LMI objects.

% We have to save the persistent variables in the SDPVAR class
obj.savedata.internalsdpvarstate = yalmip('getinternalsdpvarstate');
obj.savedata.internalsetstate = yalmip('getinternalsetstate');
obj.savedata.version = yalmip('version');

out = obj;
