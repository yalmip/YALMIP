function X = ncvar(X)
X = struct(X);
X = rmfield(X,'originalbasis');
X = rmfield(X,'leftfactors');
X = rmfield(X,'rightfactors');
X = rmfield(X,'midfactors');
X = ncvar(struct(X));