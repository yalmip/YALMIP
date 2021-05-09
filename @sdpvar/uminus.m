function Y = uminus(Y)
%UMINUS (overloaded)

if isa(Y,'blkvar')
    Y = sdpvar(Y);
end

switch Y.typeflag
    case {0,1,2,3,4}
        Y.basis = -Y.basis;   
    case 5
        error('Cone object cannot be negated');
    otherwise
end
% Reset info about conic terms
Y.conicinfo = [0 0];
Y.extra.opname='';


