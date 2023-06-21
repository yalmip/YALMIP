function x = fix_downforce(p,x)
% ensure x1 + x2 + ... <= y
for i = 1:length(p.downForce)
    forcing = p.downForce{i}.forcing;
    if x(forcing)==0
        forced = p.downForce{i}.forced;
        x(forced)=0;
    end
end