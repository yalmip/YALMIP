function X = blkvar
% BLKVAR Constructor for block variables

% Author Johan Löfberg
% $Id: blkvar.m,v 1.3 2005-06-02 13:40:00 joloef Exp $

X.blocks = {};
X = class(X,'blkvar',sdpvar(1));
	