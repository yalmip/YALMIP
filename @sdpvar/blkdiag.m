function y = blkdiag(varargin)
%BLKDIAG (overloaded)

if nargin<2
    y=varargin{1};
    return
end

% Get dimensions
n = zeros(length(varargin),1);
m = zeros(length(varargin),1);
isasdpvar = zeros(length(varargin),1);
%Symmetric = zeros(nargin,1);
for i = 1:length(varargin)
    if isa(varargin{i},'sdpvar')
        isasdpvar(i) = 1;
        n(i)=varargin{i}.dim(1);
        m(i)=varargin{i}.dim(2);
    else
        [n(i) m(i)] = size(varargin{i});
    end
end

% Find all free variables used
all_lmi_variables = [];
for i = 1:length(varargin)
    all_lmi_variables = [all_lmi_variables getvariables(varargin{i})];
end
all_lmi_variables = unique(all_lmi_variables);

% Create an SDPVAR
y=sdpvar(1,1,'rect',all_lmi_variables,[]);


% Some indexation tricks
msums = cumsum([0; m]);
nsums = cumsum([0; n]);
summ=sum(m);
sumn=sum(n);
indextable = reshape(1:sumn*summ,sumn,summ);

is = [];
js = [];
ss = [];
for j = 1:length(varargin)
    nnindex = indextable(1+nsums(j):nsums(j+1),1+msums(j):msums(j+1));
    if isasdpvar(j)
        try
            this_uses = find(ismembc(all_lmi_variables,varargin{j}.lmi_variables));            
        catch
            % Octave fix...
            this_uses = find(ismember(all_lmi_variables,varargin{j}.lmi_variables));
        end
        mindex = [1 this_uses+1];

        [a,b,d] = find(varargin{j}.basis.');
        is = [is(:);reshape(mindex(a),[],1)];
        js = [js(:);reshape(nnindex(b),[],1)];
        ss = [ss(:);d(:)];
    else
        [a,b,d] = find( varargin{j}(:).');
        is = [is;ones(length(a),1)];
        js = [js;reshape(nnindex(b),[],1)];
        ss = [ss(:);d(:)];       
    end
end

y.basis = sparse(js,is,ss,sum(m)*sum(n),1+length(all_lmi_variables));
y.dim(1) = sumn;
y.dim(2) = summ;
% Reset info about conic terms
y.conicinfo = [0 0];


