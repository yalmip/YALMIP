function AConstraint = assignschur(AConstraint,thecompiler,varargin)
AConstraint.clauses{1}.schurfun  = thecompiler;
AConstraint.clauses{1}.schurdata = varargin;