function pnew = degreduce(p,d)
% DEGREDUCE Remove higher order terms

% Author Johan Löfberg
% $Id: degreduce.m,v 1.3 2004-10-04 07:51:13 johanl Exp $

if islinear(p)
    pnew = p;
else
    [sqrList,CompressedList] = yalmip('nonlinearvariables');

    base = getbase(p);
    vars = getvariables(p);

    for i = 1:length(vars)
        v = vars(i);
        if ismember(v,CompressedList(:,1))
            j = find(CompressedList(:,1)==v);
            if sum(any(CompressedList(j,2:end),1))>d
                base(i+1)=0;
            end
        end
    end
    pnew = clean(sdpvar(1,1,[],vars,base));

end