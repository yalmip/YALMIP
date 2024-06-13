function p = removeCuts(p)

% Remove EQ cuts
p.F_struc(p.KCut.f,:)=[];
p.K.f = p.K.f - length(p.KCut.f);
p.KCut.f = [];

% Remove LP cuts
p.F_struc(p.K.f+p.KCut.l,:)=[];
p.K.l = p.K.l - length(p.KCut.l);
p.KCut.l = [];

% Remove SDP cuts
p = removeSDPcuts(p);    