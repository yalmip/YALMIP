function X = not(X)
% Internal class for constraint list

superiorto('sdpvar');
superiorto('double');

% Try to evaluate

if length(X.List)>3
    error('Negation can only be applied to BINARY relations')
else
    switch X.List{2}
        case '<'
            X.List{2} = '>';
            X.Evaluated{1} = -X.Evaluated{1};
        case '>'
            X.List{2} = '<';
            X.Evaluated{1} = -X.Evaluated{1};
        case '=='
            % A simple inequality has to be converted to an expression which
            % most likely will be a mixed-integer model....
            X = X.List{1} ~= X.List{3};       
        otherwise
            error('Negation cannot be applied on this operator')
    end
end