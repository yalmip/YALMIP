function info = sizestring(d)

info = num2str(d(1));
for i = 2:length(d)
    info = [info 'x' num2str(d(i))];
end