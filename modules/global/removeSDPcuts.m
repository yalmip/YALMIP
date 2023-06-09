function p = removeSDPcuts(p)

if length(p.KCut.s)>0
    starts = p.K.f+p.K.l + [1 1+cumsum((p.K.s).^2)];
    remove_these = [];
    for i = 1:length(p.KCut.s)
        j = p.KCut.s(i);
        remove_these = [remove_these;(starts(j):starts(j+1)-1)'];
    end
    p.F_struc(remove_these,:)=[];
    p.K.s(p.KCut.s) = [];
end