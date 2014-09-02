function sys = assignschur(AConstraint,thecompiler,varargin)

sys = sdpvar(1);
sys.typeflag = 30;
Schurfun_defined = 0;
Var_defined = 0;

for i = 1:nargin
    if isa(varargin{i},'sdpvar')
        sys.extra.variables{i} = varargin{i};
        Var_defined = 1;
    elseif isa(varargin{i},'double')
        sys.extra.data{i} = varargin{i};
    elseif isa(varargin{i},'char')
        sys.extra.Schurfun = varargin{i};
        Schurfun_defined = 1;
    end
end
sys.extra.negated = 0;

if ~Schurfun_defined
    error('A Schur-function must be defined.')
end

sys = constraint(sys,'>=',0);