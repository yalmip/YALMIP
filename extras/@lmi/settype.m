function  x = settype(F)

F = flatten(F);
nlmi = length(F.LMIid);
lmiinfo{1} = 'sdp';
lmiinfo{2} = 'elementwise';
lmiinfo{3} = 'equality';
lmiinfo{4} = 'socc';
lmiinfo{5} = 'rsocc';
lmiinfo{7} = 'integer';
lmiinfo{8} = 'binary';
lmiinfo{9} = 'kyp';
lmiinfo{10} = 'eig';
lmiinfo{11} = 'sos';
if (nlmi == 0)
    x = 'empty';
elseif nlmi == 1;
    x = lmiinfo{F.clauses{1}.type};
else
    x = lmiinfo{F.clauses{1}.type};
    for i = 1:length(F.clauses)
        if F.clauses{1}.type ~=F.clauses{i}.type
            x = 'multiple';
            break
        end
    end
end