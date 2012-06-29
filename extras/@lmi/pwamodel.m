function [A,b] = pwamodel(F,x)

[dummy,aux,s,model] = export(F,[]);

A = -model.F_struc(:,2:end);
b = model.F_struc(:,1);

keep_these = find(ismember(aux.used_variables,getvariables(x)));

A = [A(:,keep_these) A(:,setdiff(1:size(A,2),keep_these))];
[A,b] = fourier_motzkin(A,b,length(keep_these));

function [A,b] = fourier_motzkin(A,b,m)

while size(A,2)>m
    [A,b] = fourier_motzkin_1(A,b,m);
    [aux,i] = unique([A b],'rows');
    A = A(i,:);
    b = b(i);
end

function [Aout,bout] = fourier_motzkin_1(A,b,m)

for i = m+1:size(A,2)
    less = find(A(:,i)>0);
    larger = find(A(:,i)<0);
    t(i-m) =length(less)*length(larger);
end

[minn,remove] = min(t);
remove = remove+m;
keep = setdiff(1:size(A,2),remove);

less = find(A(:,remove)>0);
larger = find(A(:,remove)<0);
notinvolved = find(A(:,remove)==0);

Aout = A(notinvolved,keep);
bout = b(notinvolved);

for i = 1:length(less)
    for j = 1:length(larger)
        Aout = [Aout;A(less(i),keep)/abs(A(less(i),remove)) + A(larger(j),keep)/abs(A(larger(j),remove))];
        bout = [bout;b(less(i))+b(larger(j))];
    end
end