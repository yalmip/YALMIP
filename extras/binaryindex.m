function [z,F] = binaryindex(n,m,Fin)

% Number of bits needed
k=ceil(log2(m));

w=2^k;
b=dec2bin(0:w-1);
for i = 1:m
    pattern{i} = str2num(strrep(strrep(b(i,:),'1',' 1 '),'0',' 0 '));
end

z=sdpvar(n,m,'full');
delta = binvar(n,k,'full');
F=([]);
M=m;
for iz=1:n
    for jz=1:m
        t=0;
        for i = 1:length(pattern{jz})
            if pattern{jz}(i)
                F=F+( -M*delta(iz,i) <= z(iz,jz) <= M*delta(iz,i));
                t=t+1-delta(iz,i);
            else
                F=F+( -M*(1-delta(iz,i)) <=z(iz,jz)<= M*(1-delta(iz,i)));
                t=t+delta(iz,i);
            end
        end
        F=F+(-M*t <= z(iz,jz)-1 <= M*t);       
    end
end
F=F+(m-1 >= delta*2.^((k-1:-1:0)'));

if nargin == 3
    F= F+Fin;
end