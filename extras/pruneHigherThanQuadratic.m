function [c,p]= pruneHigherThanQuadratic(c,p)
% Used by sdpvar/expectation
% First prune non-used
keep = true(1,size(c,2));
for i = 1:length(keep)
    ci = c(:,i);
    if isa(ci,'double') && nnz(ci)==0
        keep(i) = false;
    end
end
c = c(:,keep);
p = p(keep);
% % Now prune the high-degrees
% keep = true(1,length(p));
% for i = 1:length(p)
%     if degree(p(i)) > 2
%         keep(i) = false;
%     end
% end
% if ~all(keep)
%     p = p(keep);
%     if isa(c,'cell')
%         
%     else
%         c = c(:,keep);
%     end
% end