function [isconvexorconcave,A,B] = isconvexmeshdata(xi,yi,zi)

[n,m] = size(xi);

% Go through all triangles
isconvex = 1;
isconcave = 1;
A = [];
B = [];
for i = 1:n-1
    for j = 1:m-1
        x1 = [xi(i,j);yi(i,j);zi(i,j)];
        x2 = [xi(i,j+1);yi(i,j+1);zi(i,j+1)];
        x3 = [xi(i+1,j);yi(i+1,j);zi(i+1,j)];
        
        temp = null([[x1';x2';x3'] ones(3,1)]);
        a = temp(1:3);
        b = temp(4);
        if a(3) < 0
            a = -a;
            b = -b;
        end
        isconvex = isconvex && all(a'*[xi(:)';yi(:)';zi(:)']+b >= -1e-10);
        isconcave = isconcave && all(a'*[xi(:)';yi(:)';zi(:)']+b <= 1e-10);
        A = [A a];
        B = [B b];
        
        x3 = [xi(i+1,j+1);yi(i+1,j+1);zi(i+1,j+1)];
        x2 = [xi(i,j+1);yi(i,j+1);zi(i,j+1)];
        x1 = [xi(i+1,j);yi(i+1,j);zi(i+1,j)];
        
        temp = null([[x1';x2';x3'] ones(3,1)]);
        a = temp(1:3);
        b = temp(4);
        if   a(3) < 0
            a = -a;
            b = -b;
        end     
        isconvex = isconvex && all(a'*[xi(:)';yi(:)';zi(:)']+b >= -1e-10);
        isconcave = isconcave && all(a'*[xi(:)';yi(:)';zi(:)']+b <= 1e-10);
        A = [A a];
        B = [B b];
        
        if ~isconvex && ~isconcave;
            break
        end
    end
end
% 0 : Neither, 1 Convex -1 Concave
isconvexorconcave = isconvex-isconcave;