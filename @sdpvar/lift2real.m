function X = lift2real(X)

% Given (ar+i*ai)*X*(br+i*bi) 
% Silly code, since it is Hermitian, we know ar = br', ai=-bi'
reF=real(X);
imF=imag(X);
newleftfactors = {};
newmidfactors = {};
newrightfactors = {};
for i = 1:length(X.midfactors)
    ar = real(X.leftfactors{i});
    ai = imag(X.leftfactors{i});
    br = real(X.rightfactors{i});
    bi = imag(X.rightfactors{i});
    newleftfactors{end+1} = [ar;ai];
    newmidfactors{end+1} = X.midfactors{i};
    newrightfactors{end+1} = [br -bi];
    newleftfactors{end+1} = [-ai;ar];
    newmidfactors{end+1} = X.midfactors{i};
    newrightfactors{end+1} = [bi br];
end
X.leftfactors = [];
X.midfactors = [];
X.rightfactors = [];
X = [reF -imF;imF reF];
%X = kron(eye(2),reF) + kron([0 -1;1 0],imF);

X.leftfactors = newleftfactors;
X.midfactors = newmidfactors;
X.rightfactors = newrightfactors;