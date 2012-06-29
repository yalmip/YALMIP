function sys = lazybasis(n,m,i,j,s)

%sys.basis = basis;
sys.n = n;
sys.m = m;
sys.iX = i;
sys.jX = j;
sys.sX = s;
sys = class(sys,'lazybasis');
