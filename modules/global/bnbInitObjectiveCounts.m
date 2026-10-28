function [s,info] = bnbInitObjectiveCounts(p,upper)
% Keep the count tuple behind each objective interval, with bounded work.
[s,info]=describeObjectiveCounts(p);
timer=tic;
if isempty(s) || ~isfinite(upper), s=[]; return; end
ng=numel(s.weights);
counts=accumarray(s.groups,s.roothi-s.rootlo,[ng 1]);
cap=upper-s.base+s.radius;
values=0; tuples=zeros(1,ng);
for k=1:ng
    maxCount=min(counts(k),max(0,floor((cap+s.radius)/s.weights(k))));
    if numel(values)*(maxCount+1)>65536
        s=[]; info.reason='tuple work budget'; info.seconds=toc(timer); return
    end
    next=values(:)+(0:maxCount)*s.weights(k);
    nextTuples=repmat(tuples,maxCount+1,1);
    nextTuples(:,k)=reshape(repmat(0:maxCount,numel(values),1),[],1);
    keep=next(:)<=cap+s.radius;
    values=next(keep); values=values(:); tuples=nextTuples(keep,:);
    if numel(values)>4096
        s=[]; info.reason='tuple size budget'; info.seconds=toc(timer); return
    end
end
[values,order]=sort(values);
s.values=s.base+values; s.tuples=tuples(order,:); s.cap=upper;
s.allVariables=find(p.c); s.allRootLo=ceil(p.lb(s.allVariables));
s.allRootHi=floor(p.ub(s.allVariables));
s.countRows=sparse(s.groups,1+s.variables,1-2*s.negative,ng,numel(p.c)+1);
offset=-s.rootlo; offset(s.negative)=s.roothi(s.negative);
s.countRows(:,1)=accumarray(s.groups,offset,[ng 1]);
info.reason='active'; info.entries=numel(values); info.seconds=toc(timer);
end

function [s,info] = describeObjectiveCounts(p)
% Describe a discrete objective using a few coefficient classes.
s = [];
info = struct('reason','unsupported','entries',0,'classes',0,'seconds',0);
timer = tic;
if nnz(p.Q) || ~isempty(p.evalMap) || any(p.variabletype) || ...
        ~isequal(p.K.m,0) || p.IntegerCost
    return
end
used = find(p.c);
if isempty(used) || ~all(ismember(used,p.integral_variables)), return; end
c = full(p.c(used)); c = c(:);
lo = ceil(p.lb(used)); hi = floor(p.ub(used)); lo=lo(:); hi=hi(:);
if any(~isfinite([lo;hi;c])) || any(lo>hi) || ...
        any(abs([lo;hi])>1e6) || numel(c)>2000
    return
end
basePoint = lo; basePoint(c<0) = hi(c<0);
base = p.f+c'*basePoint;
width = hi-lo;
scale = 1+abs(p.f)+sum(abs(c).*max(abs(lo),abs(hi)));
rounding = 128*eps*max(1,numel(c))*scale;
if ~isfinite(scale), return; end
active = find(width>0);
[weights,order] = sort(abs(c(active))); active=active(order);
representatives = []; groups = zeros(numel(active),1); err=0;
for i=1:numel(active)
    if isempty(representatives) || weights(i)-representatives(end)>1e-10*max(1,weights(i))
        representatives(end+1,1)=weights(i); %#ok<AGROW>
    end
    groups(i)=numel(representatives);
    err=err+width(active(i))*abs(weights(i)-representatives(end));
end
info.classes=numel(representatives);
if info.classes>8 || isempty(active)
    info.reason='coefficient classes'; return
end
s=struct('variables',used(active),'negative',c(active)<0, ...
    'rootlo',lo(active),'roothi',hi(active),'weights',representatives, ...
    'groups',groups,'base',base,'radius',err+rounding);
info.reason='active'; info.seconds=toc(timer);
end
