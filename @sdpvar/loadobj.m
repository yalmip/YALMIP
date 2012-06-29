function out = loadobj(obj)
%LOADOBJ (overloaded)

% Restore persistent variables in SDPVAR class
if ~isempty(obj.savedata)
    try
        % Some backwards compatability        
        if ~isfield(obj.savedata,'conicinfo')
            obj.conicinfo = [0 0];
        end
        if isfield(obj.savedata,'optsolution')
            obj.savedata.internalsdpvarstate.optSolution = obj.savedata.optsolution;
        end
        if isfield(obj.savedata.internalsdpvarstate,'optSolution')
            if ~isfield(obj.savedata.internalsdpvarstate.optSolution,'values')
                obj.savedata.internalsdpvarstate.optSolution.values = [];
            end
        end
        if ~isfield(obj.savedata.internalsdpvarstate,'evalVariables')           
            obj.savedata.internalsdpvarstate.evalVariables = [];            
        end
        if ~isfield(obj.savedata.internalsdpvarstate,'logicVariables')           
            obj.savedata.internalsdpvarstate.logicVariables = [];            
        end        
        if ~isfield(obj.savedata.internalsdpvarstate,'complexpair')           
            obj.savedata.internalsdpvarstate.complexpair = [];            
        end
        
        
        yalmip('setinternalsdpvarstate',obj.savedata.internalsdpvarstate);
    catch
        error('Data probably saved in old YALMIP version. Cannot load this...');
    end
end
if ~isfield(obj.savedata,'conicinfo')
   obj.conicinfo = [0 0];
end
if ~isa(obj,'sdpvar'),   
    out = class(obj,'sdpvar');
else
    out = obj;
end
out.savedata = [];