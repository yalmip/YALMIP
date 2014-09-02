function AConstraint = assignschur(AConstraint,thecompiler,varargin)
AConstraint = flatten(AConstraint);
AConstraint.clauses{1}.schurfun  = thecompiler;
AConstraint.clauses{1}.schurdata = varargin;