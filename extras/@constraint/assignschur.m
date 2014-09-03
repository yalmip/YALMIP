function sys = assignschur(AConstraint,thecompiler,varargin)

sys = assignschur(lmi(AConstraint),thecompiler,varargin{:})
