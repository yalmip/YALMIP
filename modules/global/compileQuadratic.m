function [Q,c] = compileQuadratic(c,p,onlysquares)
Q = spalloc(length(c),length(c),0);
%c = p.c;
for i = 1:size(p.bilinears,1)
    quadratic = (p.bilinears(i,2)==p.bilinears(i,3));
    if onlysquares == 3 && quadratic
        if c(p.bilinears(i,1))>0
            Q(p.bilinears(i,2),p.bilinears(i,3)) = c(p.bilinears(i,1));            
            c(p.bilinears(i,1)) = 0;
        end
    elseif c(p.bilinears(i,1)) && (~onlysquares || quadratic) && ~(onlysquares == 3)
        %if ~quadratic || (quadratic && c(p.bilinears(i,1))>0)
            Q(p.bilinears(i,2),p.bilinears(i,3)) = c(p.bilinears(i,1))/2;
            Q(p.bilinears(i,3),p.bilinears(i,2)) = Q(p.bilinears(i,3),p.bilinears(i,2))+c(p.bilinears(i,1))/2;
            c(p.bilinears(i,1)) = 0;
       % end
    end
end