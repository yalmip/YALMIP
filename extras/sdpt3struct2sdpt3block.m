function [C,A,b,blk] = sdpt3struct2sdpt3block(F_struc,c,K)
%SDPT3STRUCT2SDPT3BLOCK Internal function to convert data to SDPT3 format

% Author Johan Löfberg
% $Id: sdpt3struct2sdpt3block.m,v 1.2 2004-07-02 08:17:31 johanl Exp $


nvars = size(F_struc,2)-1;
block  = 1;
top = 1;
block = 1;
blksz = 100;

if K.l>0
	blk{block,1} = 'l';
	blk{block,2} = K.l;
	C{block,1} = sparse(F_struc(top:top+K.l-1,1));
	for var_index = 1:nvars
		A{block,var_index} =  -sparse(F_struc(top:top+K.l-1,var_index+1));
	end
	block = block+1;
	top = top+sum(K.l);
end

if K.q>0
	blk{block,1} = 'q';
	blk{block,2} = K.q;
	C{block,1} = sparse(F_struc(top:top+sum(K.q)-1,1));
 	for var_index = 1:nvars
		A{block,var_index} =  -sparse(F_struc(top:top+sum(K.q)-1,var_index+1));
	end   
	block = block+1;
	top = top+sum(K.q);
end

if K.s>0
	constraints = 1;
	while constraints<=length(K.s)
		n = K.s(constraints);
		Cvec = F_struc(top:top+n^2-1,1);
		C{block,1} = reshape(Cvec,n,n);
		for var_index = 1:nvars
			Avec = -F_struc(top:top+n^2-1,var_index+1);
			A{block,var_index} = reshape(Avec,n,n);
		end
		blk{block,1} = 's';
		blk{block,2} = n;
		top = top+n^2;
		constraints = constraints+1;
		sum_n = n;
		while (sum_n<blksz) & (constraints<=length(K.s))
			n = K.s(constraints);
			Cvec = F_struc(top:top+n^2-1,1);
			C{block,1} = blkdiag(C{block,1},reshape(Cvec,n,n));
            [n1,m1] = size(A{block,var_index});
            [n2,m2] = size(reshape(-F_struc(top:top+n^2-1,1+1),n,n));
            Z = spalloc(n1,m2,0);
			for var_index = 1:nvars
				Avec = -F_struc(top:top+n^2-1,var_index+1);
				A{block,var_index} = [A{block,var_index} Z;Z' reshape(Avec,n,n)];
			end
			blk{block,2} = [blk{block,2} n];
			top = top+n^2;
			sum_n = sum_n+n;
			constraints = constraints+1;
		end
		block = block+1;
	end
end

% And we solve dual...
b = -c(:);

% blkdiag with 2 sparse blocks
function y = blkdiag(x1,x2)
[n1,m1] = size(x1);
[n2,m2] = size(x2);
Z = spalloc(n1,m2,0);
y = [x1 Z;Z' x2];