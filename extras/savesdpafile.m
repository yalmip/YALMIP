function solution = savesdpafile(varargin)
%SAVESDPAFILE Saves a problem definition in the SDPA format
%
%    SAVESDPAFILE(F,h,'filename')    Saves the problem min(h(x)), F(x)>0 to the file filename
%    SAVESDPAFILE(F,h)               A "Save As" - box will be opened                               
%
% Note the the SDPA format does not support SOCPs or equalities.
% Equalities will be eliminated using double-sided inequalities. If the
% model contains SOCP constraints the command will exit.

F = varargin{1};
h = varargin{2};
nvars = yalmip('nvars');

if isa(F,'constraint')
    F = lmi(F);
end

if any(is(F,'socp'))
	error('savesdpafile does not support SOCPs (not supported in the SDPA format');
end

% Expand nonlinear operators
[F,failure,cause] = expandmodel(F,h,sdpsettings);
if failure % Convexity propgation failed
    interfacedata = [];
    recoverdata = [];
    solver = '';
    diagnostic.solvertime = 0;
    diagnostic.problem = 14;
    diagnostic.info = yalmiperror(14,cause);
    return
end

% Convert equalities to inequalities
feq =find(is(F,'equality'));
if ~isempty(feq)
    f = sdpvar(F(feq));
    F(feq)=[];
    F = [F, f >= 0, f <= 0];
end

% Get the SP format
[F_struc,K]  = lmi2sedumistruct(F);

% Convert the objective
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

if K.f>0
	% Extract the inequalities
	A_equ = F_struc(1:K.f,2:end);
	b_equ = -F_struc(1:K.f,1);
	[Q,R] = qr(A_equ');
	n = max(find(any(R')));
	Q1 = Q(:,1:n);
	Q2 = Q(:,n+1:end); 
	R = R(1:n,:);
	x_equ = Q1*(R'\b_equ);  
	% Exit if no consistent solution exist
	if (norm(A_equ*x_equ-b_equ)>sqrt(eps))
		error('Linear constraints inconsistent.');
		return
	end
	% So we dont need these rows anymore
	F_struc = F_struc(K.f+1:end,:);
	K.f = 0;
	% OK, we found a new basis
	H = Q2;
	% objective in new basis
	c = H'*c;
	% LMI in new basis
	F_struc = [F_struc*[1;x_equ] F_struc(:,2:end)*H];
end

% Is a filename supplied
if nargin<3
	[filename, pathname] = uiputfile('*.dat-s', 'Save SDPA sparse format file');
	if isa(filename,'double')
		return % User canceled
	else
		% Did the user change the extension
		if isempty(strfind(filename,'.'))
			filename = [pathname filename '.dat-s'];
		else
			filename = [pathname filename];
		end	
	end
else
	filename = varargin{3};
end

% Save to file
integer_variables = find(ismember(used_variables,yalmip('intvariables')));
createsdplibfile(F_struc, K, c, filename,integer_variables);
