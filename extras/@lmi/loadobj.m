function out = loadobj(obj)
%loadOBJ load filter for LMI objects.

% Restore persistent variables in SDPVAR class
if ~isempty(obj.savedata)
    try
        % Some backwards compatability
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
        yalmip('setinternalsetstate',obj.savedata.internalsetstate);
    catch
        error('Data probably saved in old YALMIP version. Cannot load this...');
    end
end
if ~isa(obj,'lmi'),    
    if ~isfield(obj,'LMIid')
        LMIid = [];
        obj = rmfield(obj,'savedata');
        for i = 1:length(obj.clauses)
            LMIid = [LMIid obj.clauses{i}{1}.LMIid];
            X.data=obj.clauses{i}{1}.data;
            X.type = obj.clauses{i}{1}.type;
            X.symbolic = obj.clauses{i}{1}.symbolic;
            X.handle = obj.clauses{i}{1}.handle;
            X.strict = obj.clauses{i}{1}.strict;
            X.cut = obj.clauses{i}{1}.cut;
            obj.clauses{i} = X;
        end
        obj.LMIid = LMIid;
    end    
    obj.savedata = [];
    out = class(obj,'lmi');
else
    out = obj;
end
out.savedata = [];