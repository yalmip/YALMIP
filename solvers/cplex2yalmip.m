function [H,C,A,B,LB,UB,QC,VARTYPE,INDEQ,PARAM,OPTIONS] = cplex2yalmip(interfacedata)

options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
Q       = interfacedata.Q;
K       = interfacedata.K;
x0      = interfacedata.x0;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
UB      = interfacedata.ub;
LB      = interfacedata.lb;

SENSE = 1;     
C = full(c);   
if K.l+K.f+K.q == 0
    A = zeros(1,length(c));A(1)=1;
    B = 1e6;
else
    A =-(F_struc(1:K.f+K.l,2:end)); 
    B = full(F_struc(1:K.f+K.l,1));            
end

INDEQ = [];
if K.f>0
    INDEQ(1:K.f) = 1:K.f;
end

if any(K.q)
    top = K.f+K.l + 1;
    for i = 1:length(K.q)
        % [cx+d;Ax+b]   |Ax+b|<cx+d, originally a QCQP
        m = K.q(i);
        ci = F_struc(top,2:end)';
        di = F_struc(top,1);
        Ai = F_struc(top+1:top+m-1,2:end);
        bi = F_struc(top+1:top+m-1,1);        
        QC(i).Q = full(Ai'*Ai - ci*ci');
        QC(i).r = full(di'*di - bi'*bi);   
        QC(i).L = full(2*bi'*Ai - 2*di*ci');
        top = top+m;
    end
else
    QC = [];
end

VARTYPE = repmat('C',size(A,2),1);
VARTYPE(integer_variables)='I'; % Integer variables
VARTYPE(binary_variables) ='B'; % Binary variables

if nnz(Q)==0
    H = [];
else
    H = full(2*Q);
end

try
    PARAM = options.cplex.param;
catch
    PARAM = [];
end
OPTIONS.verbose = options.verbose;
try
OPTIONS.logfile = options.cplex.logfile;
catch
    OPTIONS.logfile = '';
end
if ~isempty(x0)
    OPTIONS.x0 = [(1:length(x0))' x0(:)];
end