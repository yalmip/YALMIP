function X = blkvar
% BLKVAR Constructor for block variables

X.blocks = {};
X = class(X,'blkvar',sdpvar(1));
	