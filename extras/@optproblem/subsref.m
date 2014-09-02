function P = subsref(P,String)
%subsref           Overloaded indexing

if length(String)>1
    R = subsref(P,String(1));
    P = subsref(R,String(2:end));
    return
end
switch String.type
    case '.'
        switch String.subs
            case 'Constraints'
                P = P.Constraints;
            case 'Objective'
                P = P.Objective;
            case 'Options'
                P = P.Options;
            case 'maximize'
                maximize(P);
            case 'minimize'
                minimize(P);
            otherwise
                error('Indexing type not supported');
        end
    case '()'
        if isa(String.subs{1},'struct')
            P.Options = String.subs{1};
        elseif isa(String.subs{1},'lmi') | isa(String.subs{1},'constraint')
            P = [P, String.subs{1}];
        else
            error('Indexing type not supported');
        end
    otherwise
        error('Indexing type not supported');
end



