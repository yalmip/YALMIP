function model = yalmip2scs(interfacedata)

%% From internal
data.A = -interfacedata.F_struc(:,2:end);
data.b = full(interfacedata.F_struc(:,1));
data.c =  interfacedata.c;
cones = [];
cones.f = interfacedata.K.f;
cones.l = interfacedata.K.l;
cones.q = interfacedata.K.q;
cones.s = interfacedata.K.s;
cones.ep =interfacedata.K.e;
param = interfacedata.options.scs;
param.verbose = interfacedata.options.verbose;

%% Extract lower diagonal form for new SCS format
if ~isempty(cones.s) && any(cones.s)
    sdpA = data.A(1+cones.l + cones.f+sum(cones.q):end,:);
    sdpb = data.b(1+cones.l + cones.f+sum(cones.q):end,:);
    expA = data.A(end-3*cones.ep+1:end,:);
    expb = data.b(end-3*cones.ep+1:end,:);
    data.A = data.A(1:cones.l + cones.f+sum(cones.q),:);    
    data.b = data.b(1:cones.l + cones.f+sum(cones.q),:);
    top = 1;
    for i = 1:length(cones.s)
        A = sdpA(top:top + cones.s(i)^2-1,:);
        b = sdpb(top:top + cones.s(i)^2-1,:);
        n = cones.s(i);
        ind = find(speye(n));
        b(ind) = b(ind)/sqrt(2);
        A(ind,:) = A(ind,:)/sqrt(2);
        ind = find(tril(ones(n)));
        A = A(ind,:);
        b = b(ind);
        data.A = [data.A;A];
        data.b = [data.b;b];
        top = top  + cones.s(i)^2;
    end
    data.A = [data.A;expA];
    data.b = [data.b;expb];
end

%% Collect in one structure
model.data = data;
model.cones = cones;
model.param = param;
