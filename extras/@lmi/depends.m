function LinearVariables = depends(varargin)

if nargin > 1
    LinearVariables = [];
    for i = 1:nargin
        LinearVariables_i = depends(varargin{i});
        LinearVariables = [LinearVariables;LinearVariables_i(:)];
    end
    LinearVariables = unique(LinearVariables);
    return
else
    F = varargin{1};
end
        
F = flatten(F);

% Get all used variables in this LMI object
used = recursivedepends(F.clauses);

% Now, determine the involved linear variables
[mt,variabletype] = yalmip('monomtable');
if any(variabletype)%1%~isempty(nlv) & any(ismembc(used,nlv(1,:)))
    LinearVariables = find(any(mt(used,:),1));    
    LinearVariables = LinearVariables(:)';
else
    LinearVariables = used;
end

function used = recursivedepends(clauses)

if length(clauses) > 2
    mid = floor(length(clauses)/2);
    used1 = recursivedepends({clauses{1:mid}});
    used2 = recursivedepends({clauses{mid+1:end}});
    used = uniquestripped([used1 used2]);
else
    used = [];
    for i = 1:length(clauses)        
        if clauses{i}.type == 56
            % Meta constraint such as implies. This object is just holding
            % the data involved 
            Fivar = [];
            for j = 1:length(clauses{i}.data)
                Fivar = [Fivar getvariables(clauses{i}.data{j})];
            end
        else
            Fivar = getvariables(clauses{i}.data);
        end
        used = uniquestripped([used Fivar(:)']);
    end
end