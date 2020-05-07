function tests = test_ndsdpvar_sudoku
tests = functiontests(localfunctions);

function test1(dummy)

S = [0,0,1,9,0,0,0,0,8;6,0,0,0,8,5,0,3,0;0,0,7,0,6,0,1,0,0;...
     0,3,4,0,9,0,0,0,0;0,0,0,5,0,4,0,0,0;0,0,0,0,1,0,4,2,0;...
     0,0,5,0,7,0,9,0,0;0,1,0,8,4,0,0,0,7;7,0,0,0,0,9,2,0,0];

p = 3;
A = ndsdpvar(p^2,p^2,p^2,'full');
F = (binary(A));
F = F + (sum(A,1) == 1);
F = F + (sum(A,2) == 1);
F = F + (sum(A,3) == 1);

for m = 1:p
    for n = 1:p
        for k = 1:9
            s = sum(sum(A((m-1)*p+(1:p),(n-1)*p+(1:p),k)));  
            F = F + (s == 1);
        end
    end
end

for i = 1:p^2 
    for j = 1:p^2 
        if S(i,j)
            F = F + (A(i,j,S(i,j)) == 1);
        end
    end
end

sol = optimize(F);

Z = 0;
for i = 1:p^2
      Z = Z  + i*value(A(:,:,i));
end
assert(sol.problem == 0);
assert(norm(sort(Z(:,1)) - (1:p^2)') <= 1e-4);

