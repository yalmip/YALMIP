function sys = logdet(P)

% Create an object
if isa(P,'sdpvar')
    if is(P,'hermitian')
        superiorto('double')
        superiorto('sdpvar')
        sys.P  = {P};
        sys.cx = [];
        sys.gain = 1;
        sys = class(sys,'logdet');
    else
        error('logdet can only be applied to Hermitian SDPVAR objects')
    end
elseif isa(P,'double')
    if isessentiallyhermitian(P)
        sys = sum(log(abs(real(eig(P)))));
    else
        error('logdet can only be applied to Hermitian objects')
    end
end


