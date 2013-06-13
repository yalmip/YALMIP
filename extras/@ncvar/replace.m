function Z = replace(X,Y,W)
%REPLACE Substitutes variables
%
%Z = REPLACE(Y,X,W)  Replaces any occurence of variable Y in X with W

Z = ncvar_replace(X,Y,W);
