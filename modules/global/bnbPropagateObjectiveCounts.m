function [p,cuts,info] = bnbPropagateObjectiveCounts(s,p,lower,upper)
% Propagate bounded integer objective counts within this node only.
% Each tuple represents a cost interval, including coefficient error. Keep
% every tuple whose interval overlaps [lower,upper] and whose group counts
% fit the node domains. Their count bounds preserve all improving integer
% solutions. Returned cuts are local: they must not enter the global pool.
cuts=[];
info=struct('changed',0,'fixed',0,'pruned',0,'groups',0,'survivors',0);
if isempty(s) || ~p.feasible || isnan(lower) || ~isfinite(upper), return; end
margin=1e-7*max([1 abs(upper) abs(lower(isfinite(lower)))]);
if upper>s.cap+margin || ...
        any(ceil(p.lb(s.allVariables))<s.allRootLo) || ...
        any(floor(p.ub(s.allVariables))>s.allRootHi)
    return
end
lo=ceil(p.lb(s.variables)); hi=floor(p.ub(s.variables)); lo=lo(:); hi=hi(:);
low=lo-s.rootlo; high=hi-s.rootlo;
low(s.negative)=s.roothi(s.negative)-hi(s.negative);
high(s.negative)=s.roothi(s.negative)-lo(s.negative);
if any(low>high), p.feasible=0; info.pruned=1; return; end
ng=numel(s.weights);
groupLo=accumarray(s.groups,low,[ng 1]);
groupHi=accumarray(s.groups,high,[ng 1]);
keep=s.values>=lower-margin-s.radius & s.values<=upper+margin+s.radius;
tuples=s.tuples(keep,:);
keep=all(tuples>=groupLo' & tuples<=groupHi',2);
tuples=tuples(keep,:); info.survivors=size(tuples,1);
if isempty(tuples), p.feasible=0; info.pruned=1; return; end
minCount=min(tuples,[],1)'; maxCount=max(tuples,[],1)';
info.groups=nnz(minCount>groupLo | maxCount<groupHi);
newLow=max(low,minCount(s.groups)-groupHi(s.groups)+high);
newHigh=min(high,maxCount(s.groups)-groupLo(s.groups)+low);
newLo=s.rootlo+newLow; newHi=s.rootlo+newHigh;
newLo(s.negative)=s.roothi(s.negative)-newHigh(s.negative);
newHi(s.negative)=s.roothi(s.negative)-newLow(s.negative);
info.changed=nnz(newLo>lo)+nnz(newHi<hi);
info.fixed=nnz(newLo==newHi & lo<hi);
if info.changed
    p.lb(s.variables)=max(p.lb(s.variables),newLo);
    p.ub(s.variables)=min(p.ub(s.variables),newHi);
end
% Select in the small space of group counts before assembling sparse rows.
% Coefficients are +/-1, so 0.1 is measured in integer count units. An
% unavailable reference point disables cuts, but still permits propagation.
if ~info.groups || ~isfield(p,'x0') || isempty(p.x0), return; end
reference=s.countRows*[1;p.x0(:)];
useLo=minCount>groupLo & reference-minCount < -0.1;
useHi=maxCount<groupHi & maxCount-reference < -0.1;
if ~any(useLo) && ~any(useHi), return; end
lowCut=s.countRows(useLo,:); lowCut(:,1)=lowCut(:,1)-minCount(useLo);
highCut=-s.countRows(useHi,:); highCut(:,1)=highCut(:,1)+maxCount(useHi);
cuts=[lowCut;highCut];
end
