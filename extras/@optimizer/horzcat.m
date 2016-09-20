function P = horzcat(varargin)

for i = 1:nargin
    preEvaled(i) = isa(varargin{i},'optimizer') && isempty(varargin{i}.dimin);
end

if all(preEvaled)
    keptvariables = varargin{1}.keptvariables;
    for i = 2:nargin
        if ~isequal(keptvariables,varargin{i}.keptvariables)
            error('The decision variables in the evaluated models are not the same');
        end
    end
    model = varargin{1}.model;
    for i = 2:nargin
        model.f = model.f + varargin{i}.model.f;
        model.c = model.c + varargin{i}.model.c;
        model.Q = model.Q + varargin{i}.model.Q;           
        model = mergeConstraints(model,varargin{i}.model); 
    end
    model.f = model.f/nargin;
    model.c = model.c/nargin;
    model.Q = model.Q/nargin;   
    P = varargin{1};
    P.model = model;
    return
end

if isa(A,'optimizer') & isa(B,'optimizer')
    if isequal(A.input.expression,B.input.expression) & isequal(A.output.expression,B.output.expression)
        P = optimizer([A.F,B.F],A.h+B.h,A.ops,A.input.expression,A.output.expression);
    else
        error('The optimizer objects must share parameters and decision variables');
    end
elseif isa(B,'optimizer') & (isa(A,'constraint') | isa(A,'lmi'))
    P = A;
    A = B;
    B = P;
end
P = optimizer([A.F,B],A.h,A.ops,A.input.expression,A.output.expression);


function model = mergeConstraints(model1,model2);

top1 = 1;
top2 = 1;
K = model1.K;

Data = model1.F_struc(top1:top1 + model1.K.f-1,:);
Data = [Data;model2.F_struc(top2:top2 + model2.K.f-1,:)];
top1 = top1 + model1.K.f;
top2 = top2 + model2.K.f;
K.f = model1.K.f + model2.K.f;

Data = [Data;model1.F_struc(top1:top1 + model1.K.l-1,:)];
Data = [Data;model2.F_struc(top2:top2 + model2.K.l-1,:)];
top1 = top1 + model1.K.l;
top2 = top2 + model2.K.l;
K.l = model1.K.l + model2.K.l;

Data = [Data;model1.F_struc(top1:top1 + sum(model1.K.q)-1,:)];
Data = [Data;model2.F_struc(top2:top2 + sum(model2.K.q)-1,:)];
top1 = top1 + sum(model1.K.q);
top2 = top2 + sum(model2.K.q);
K.q = [model1.K.q model2.K.q];

Data = [Data;model1.F_struc(top1:top1 + sum(model1.K.s).^2-1,:)];
Data = [Data;model2.F_struc(top2:top2 + sum(model2.K.s).^2-1,:)];
top1 = top1 + sum(model1.K.s).^2;
top2 = top2 + sum(model2.K.s).^2;
K.s = [model1.K.s model2.K.s];

model = model1;
model.F_struc = Data;
model.K = K;

