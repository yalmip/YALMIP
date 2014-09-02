function display(X)
%display           Overloaded

P = X.P;
classification = 'Logdet-term ';
[n,m] = size(P);
classification = [classification num2str(n) 'x' num2str(m)];

value = double(X);
if ~isnan(value);
    classification = [classification ' (current value: ' num2str(value) ')'];
end                  
disp(classification);