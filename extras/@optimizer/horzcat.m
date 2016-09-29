function P = horzcat(varargin)

if isempty(varargin{1}) && isa(varargin{2},'optimizer')
    P = varargin{2};
    return
end

% Are they similiar objects (same parameterization). 
% FIXME: Check will be strengthened
for i = 1:nargin
    if isa(varargin{i},'optimizer') && isequal(varargin{i}.diminOrig,varargin{1}.diminOrig)
    else
        error('The optimizer objects do not match')
    end
end

% Do they involve the same variables?
% FIXME: Allow different sets as long as relevant to merge
for i = 1:nargin
    if ~isequal(varargin{i}.model.used_variables,varargin{1}.model.used_variables)    
        error('The optimizer objects do not match')
    end
end

% Set up some hash structures to enable pruning of repeated constraints
hashvar = randn(size(varargin{1}.model.F_struc,2),1);
n = size(varargin{1}.model.F_struc,1);for i = 2:nargin;n = max(n,size(varargin{i}.model.F_struc,1));end
hashcon = randn(n,1);
for i = 1:length(varargin)
   varargin{i}.model = defineHash(varargin{i}.model,hashvar,hashcon);
end

% Merge all models
model = varargin{1}.model;
for i = 2:nargin
    model.f = model.f + varargin{i}.model.f;
    model.c = model.c + varargin{i}.model.c;
    model.Q = model.Q + varargin{i}.model.Q;
    model = mergeConstraints(model,varargin{i}.model);
end

% Use average objective.
model.f = model.f/nargin;
model.c = model.c/nargin;
model.Q = model.Q/nargin;

P = varargin{1};
P.model = model;


function model = mergeConstraints(model1,model2);
top1 = 1;
top2 = 1;
K = model1.K;
Data = [];
hash = model1.hash;

if model1.K.f > 0
    Data = model1.F_struc(top1:top1 + model1.K.f-1,:);
end
if model2.K.f > 0
    index = top2:top2 + model2.K.f-1;
    index = find(~ismember(model2.hash.f,model1.hash.f));
    if ~isempty(index)
        Data = [Data;model2.F_struc(top2:top2 + model2.K.f-1,:)];
        K.f = K.f + length(index);
        hash.f = [hash.f;model2.hash.f(index)]
    end
end
top1 = top1 + model1.K.f;
top2 = top2 + model2.K.f;

if model1.K.l > 0
    Data = [Data;model1.F_struc(top1:top1 + model1.K.l-1,:)];
    hash.l = model1.hash.l;
end
top1 = top1 + model1.K.l;
K.l = model1.K.l;

if model2.K.l > 0
    index = top2:top2 + model2.K.l-1;
    keep = find(~ismember(model2.hash.l,model1.hash.l));
    index = index(keep);
    if ~isempty(index)
        Data = [Data;model2.F_struc(index,:)];
        K.l = K.l + length(index);
        hash.l = [hash.l;model2.hash.l(keep)];
    end
end
top2 = top2 + model2.K.l;

if ~isempty(model1.K.q) && any(model1.K.q)
    Data = [Data;model1.F_struc(top1:top1 + sum(model1.K.q)-1,:)];
end
top1 = top1 + sum(model1.K.q);
K.q = model1.K.q;

if ~isempty(model2.K.q) && any(model2.K.q)        
    keep = find(~ismember(model2.hash.q,model1.hash.q));
    if any(keep)
        for i = 1:length(keep)
        if keep(i)
            Data = [Data;model2.F_struc(top2:top2 + model2.K.q(i)-1,:)];
            K.q = [K.q model2.K.q(i)];
            hash.q = [hash.q model2.hash.q(i)];
            top2 = top2 + model2.K.q(i);
        else
            top2 = top2 + model2.K.q(i);
        end
        end
    else
        top2 = top2 + sum(model2.K.q);
    end
end

if ~isempty(model1.K.s) && any(model1.K.s)
    Data = [Data;model1.F_struc(top1:top1 + sum(model1.K.s.^2)-1,:)];
end
top1 = top1 + sum(model1.K.s.^2);
K.s = model1.K.s;

if ~isempty(model2.K.s) && any(model2.K.s)        
    keep = find(~ismember(model2.hash.s,model1.hash.s));
    if any(keep)
        for i = 1:length(keep)
        if keep(i)
            Data = [Data;model2.F_struc(top2:top2 + model2.K.s(i)^2-1,:)];
            K.s = [K.s model2.K.s(i)];
            hash.s = [hash.s model2.hash.s(i)];
            top2 = top2 + model2.K.s(i)^2;
        else
            top2 = top2 + model2.K.s(i)^2;
        end
        end
    else
        top2 = top2 + sum(model2.K.s.^2);
    end
end

model = model1;
model.F_struc = Data;
model.K = K;
model.hash = hash;


function model = defineHash(model,hashvar,hashcon);

top = 1;
if model.K.f > 0
    model.hash.f = (model.F_struc(top:top+model.K.f-1,:)*hashvar);
    top = top + model.K.f;
else
    model.hash.f = [];
end

if model.K.l > 0
    model.hash.l = (model.F_struc(top:top+model.K.l-1,:)*hashvar);
    top = top + model.K.l;
else
    model.hash.l = [];
end

if ~isempty(model.K.q) && any(model.K.q)
    for i = 1:length(model.K.q)
        model.hash.q(i) = hashcon(1:model.K.q(i))'*(model.F_struc(top:top+model.K.q(i)-1,:)*hashvar);
        top = top + model.K.q(i);
    end
else
    model.hash.q = [];
end

if ~isempty(model.K.s) && any(model.K.s)
    for i = 1:length(model.K.s)
        model.hash.s(i) = hashcon(1:model.K.s(i)^2)'*(model.F_struc(top:top+model.K.s(i)^2-1,:)*hashvar);
        top = top + model.K.s(i);
    end
else
    model.hash.s = [];
end