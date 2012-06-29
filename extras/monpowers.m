function all_monoms = monpowers(n,d)

if max(d)==0
    all_monoms = [];
else
    all_monoms  = fliplr(eye(n));
    last_monoms = all_monoms;
end

for degrees = 1:1:d-1
    new_last_monoms = [];
    for variable = 1:n
        temp = last_monoms;
        temp(:,variable) = temp(:,variable)+1;
        new_last_monoms = [new_last_monoms;temp];
     %   all_monoms = [all_monoms;temp];
    end
    last_monoms = unique(new_last_monoms,'rows');
    all_monoms = [all_monoms;last_monoms];
    %all_monoms = unique(all_monoms,'rows');
end

all_monoms = [zeros(1,n);all_monoms];
all_monoms = fliplr(all_monoms);