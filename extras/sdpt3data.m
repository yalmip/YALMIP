function [C,A,b,blk] = sdpt3data(F,h)
%SDPT3DATA Internal function to convert data to SDPT3 format

if ~(isempty(F) | isa(F,'lmi'))
	help lmi
	error('First argument (F) should be an lmi object. See help text above');
end

if ~(isempty(h) | isa(h,'sdpvar'))
	help solvesdp
	error('Third argument (the objective function h) should be an sdpvar object (or empty). See help text above');
end

[ProblemString,real_data] = catsdp(F);

% This one is used a lot
nvars = sdpvar('nvars'); 

% Convert the objective
onlyfeasible = 0;
if isempty(h)
	c=zeros(nvars,1);
else  
	[n,m]=size(h);
	if ~((n==1) & (m==1))
		error('Scalar expression to minimize please.');
	else
		lmi_variables = getvariables(h);
		c  = zeros(nvars,1);
		for i=1:length(lmi_variables)
			c(lmi_variables(i))=getbasematrix(h,lmi_variables(i));
		end;
	end
end

[F_struc,K] = lmi2sedumistruct(F);

% Which sdpvar variables are actually in the problem
used_variables_LMI = find(any(F_struc(:,2:end),1));
used_variables_obj = find(any(c',1));
used_variables = uniquestripped([used_variables_LMI used_variables_obj]);

% Check for unbounded variables
unbounded_variables = setdiff(used_variables_obj,used_variables_LMI);
if ~isempty(unbounded_variables)
	% Remove unbounded variable from problem
	used_variables = setdiff(used_variables,unbounded_variables);
end

% Pick out the necessary rows
if length(used_variables)<nvars
	c = c(used_variables);
	F_struc = sparse(F_struc(:,[1 1+[used_variables]]));
end


if (K.f>0) 
	% Extract the inequalities
	A_equ = F_struc(1:K.f,2:end);
	b_equ = -F_struc(1:K.f,1);
	
	% Find feasible (turn off annoying warning on PC)
	% Using method from Nocedal-Wright book
	showprogress('Solving equalities',options.ShowProgress);
	[Q,R] = qr(A_equ');
	n = size(R,2);
	Q1 = Q(:,1:n);
	R = R(1:n,:);
	x_equ = Q1*(R'\b_equ);  
	% Exit if no consistent solution exist
	if (norm(A_equ*x_equ-b_equ)>sqrt(eps))
		error('Linear constraints inconsistent.');
		return
	end
	% We dont need the rows for equalities anymore
	F_struc = F_struc(K.f+1:end,:);
	K.f = 0;
	
	% We found a new basis
	H = Q(:,n+1:end); % New basis
	
	% objective in new basis
	c = H'*c;
	% LMI in new basis
	F_struc = [F_struc*[1;x_equ] F_struc(:,2:end)*H];	
else
	% For simpliciy we introduce a dummy coordinate change
	x_equ = 0;
	H     = 1;
end

[C,A,b,blk] = sdpt3struct2sdpt3block(F_struc,c,K);
A = svec(blk,A,ones(size(blk,1),1)); 