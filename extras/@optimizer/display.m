function display(X)
%display           Overloaded

n = 0;
for i = 1:length(X.diminOrig)
    n = n + prod(X.diminOrig{i});
end
m = 0;
for i = 1:length(X.dimoutOrig)
    m = m + prod(X.dimoutOrig{i});
end
text = ['Optimizer object with ' num2str(n) ' inputs and ' num2str(m) ' outputs.'];
text = [text ' Solver: ' upper(X.model.solver.tag)];
disp(text)
