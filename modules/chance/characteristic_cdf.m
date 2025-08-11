function varargout = characteristic_cdf(varargin)
% Internal function to supply numerical evaluation of CDF and its
% derivative based on characteristic function

switch class(varargin{1})
    
    case 'double'
        % Asked to evaluate cdf at point x
        % The function is Probability( h(x) + g'(x)*w <= 0 )
        % i.e.  Probability( z <= -h(x)) where z = g'(x)*w
        x = varargin{1}(:);
        funcs = varargin{2};
        distribution = varargin{3};
        
        phi = distribution.characteristicfunction;
        
        % Evaluate terms at x
        g = funcs.g(x);
        h = funcs.h(x);
        
        % Define char. func for linear combination
        phi = @(t) prod(phi(g(:)*t),1);        
        
        % Compute the cdf
        varargout{1} = compute_cdf_using_phi(-h, phi);
        
    case 'sdpvar'       
        varargout{1} = yalmip('define',mfilename,varargin{:});    
        
    case 'char'
        
        operator = struct('convexity','none','monotonicity','increasing','definiteness','positive','model','callback');        
        operator.range = [0 1];
        
        % Create a derivative callback 
        funcs = varargin{4};
        distribution = varargin{5};
        phi = distribution.characteristicfunction;
        dphi = distribution.characteristicfunction_derivative;
        reldphi = distribution.characteristicfunction_relativederivative;
        % Create a function which computes gradient at x
        operator.derivative = @(x)compute_dcdf_using_phi(x,funcs.h,funcs.dh,funcs.g,funcs.dg,phi,dphi,reldphi);
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};
        
    otherwise
end

function cdf = compute_cdf_using_phi(y,phi)
% Perform the inverse Fourier transform to obtain the CDF
% for Probability(z <= y) where z has characterstic function phi(t)

integrand = @(t) imag(phi(t) .* exp(-1i * t * y)./t);
cdf =  .5-integral(integrand,0, inf)/pi;


function dcdf = compute_dcdf_using_phi_finite_difference(x,h,dh,g,dg, phi,dphi)
h0 = h(x);
dh0 = dh(x);
g0 = g(x);
dg0 = dg(x);

% Silly numerical differentiation (which is what a solver does anyway)
phi0 = @(t) prod(phi(g0(:)*t),1);  
cdf0 = compute_cdf_using_phi(-h0,phi0);
eps = 1e-8;
dcdf = [];
for k = 1:length(x)
    x_ = x;x_(k) = x_(k)+eps;
    g_ = g(x_);
    h_ = h(x_);    
    phi_ = @(t) prod(phi(g_(:)*t),1);  
    cdf_ = compute_cdf_using_phi(-h_,phi_);
    dcdf = [dcdf;(cdf_-cdf0)/eps];
end

%% -------------------------------------------------------------------------------------- 
function dcdf = compute_dcdf_using_phi(x,h,dh,g,dg,phi,dphi,reldphi)
% Compute derivative of Probability(h(x)+g(x)'*w <= 0)
h0 = h(x);
dh0 = dh(x);
g0 = g(x);
dg0 = dg(x);

% fixed parameters
DEFAULT_MAXINTERVALCOUNT = 16384;
ATOL = 1e-10;
RTOL = 1e-6;
opstruct = struct('ThrowOnFail', false);

% nodes and weights for Gauss-Kronrod 7/15 method
% nodes and weights for Gauss-Kronrod 7/15 method
rule.Nodes = [ ...
    -0.9914553711208126, -0.9491079123427585, -0.8648644233597691, ...
    -0.7415311855993944, -0.5860872354676911, -0.4058451513773972, ...
    -0.2077849550078985, 0, 0.2077849550078985, ...
    0.4058451513773972, 0.5860872354676911, 0.7415311855993944, ...
    0.8648644233597691, 0.9491079123427585, 0.9914553711208126];
rule.HighWeights = [ ...
    0.02293532201052922, 0.06309209262997855, 0.1047900103222502, ...
    0.1406532597155259, 0.1690047266392679, 0.1903505780647854, ...
    0.2044329400752989, 0.2094821410847278, 0.2044329400752989, ...
    0.1903505780647854, 0.1690047266392679, 0.1406532597155259, ...
    0.1047900103222502, 0.06309209262997855, 0.02293532201052922];
rule.LowWeights = [ ...
    0, 0.1294849661688697, 0, ...
    0.2797053914892767, 0, 0.3818300505051189, ...
    0, 0.4179591836734694, 0, ...
    0.3818300505051189, 0, 0.2797053914892767, ...
    0, 0.1294849661688697, 0];
WT = rule.HighWeights;
EWT = WT - rule.LowWeights;
NODES = rule.Nodes(:);

% interval split
interval = [0,1];
pathlen = interval(end)-interval(1);

% compute phi_z(t)
phi_u = @(t) phi(g0(:)*t);
phi_z = @(t) prod(phi(g0(:)*t),1);
exp_ith = @(t) exp(1i*t*h0);

% compute f_z(-h0)

integrand1 = @(t) real(exp_ith(t) .* phi_z(t));
pdf_val = (1/pi)*integral(integrand1,0,inf);

[PHI,E,U,W,INTERVAL,pdf_value] = compute_pdf(phi_u,exp_ith,...
    interval,pathlen,WT,EWT,NODES,DEFAULT_MAXINTERVALCOUNT,ATOL,RTOL,opstruct);
pdf_value = pdf_value/pi;
DPHI = size(PHI);
for i = 1:length(g0)
    for j = 1:length(U)
        dphi_temp = dphi(g0(i)*U(j));
        DPHI(i,j) = dphi_temp(i);
    end
end

% compute dphi_z/dg
num_j = length(g0);
dphi_p1 = @(t) [];
for j = 1:num_j
    dphi_p1 = @(t) [dphi_p1(t),dphi(g0(j)*t)];
end
integrand2 = @(t) imag(exp_ith(t) .* phi_z(t) .* diag(dphi_p1(t)) ./ phi(g0(:)*t));
terms = (-1/pi)*integral(integrand2,0,inf,'ArrayValued',true);

% commpute the whole derivative
dcdf = (-pdf_val.*dh0') + terms'*dg0;

[predcdf,predcdf_errbnd] = compute_predcdf(h0,g0,dh0,dg0,exp_ith,...
    PHI,DPHI,E,W,INTERVAL,pathlen,WT,EWT,NODES,...
    DEFAULT_MAXINTERVALCOUNT,ATOL,RTOL,opstruct);
dcdf_value = (-pdf_value.*dh0') + predcdf'*dg0;

% ————————————————————————远程版本————————————————————————————
integrand_pdf = @(t) real(exp_ith(t) .* phi_z(t));

pdf_val = (1/pi)*integral(integrand_pdf,0,inf);

rr = @(t)(reldphi(g0(:)*t));
integrand_dcdf = @(t) imag(exp_ith(t) .* phi_z(t) .* rr(t));
if nnz(dg0)==0
	dcdf = (-pdf_val.*dh0');
else
    terms = (-1/pi)*integral(integrand_dcdf,0,inf,'ArrayValued',true);
    % commpute the whole derivative
    dcdf = (-pdf_val.*dh0') + terms'*dg0;
end

% dcdf = compute_dcdf_using_phi_finite_difference(x,h,dh,g,dg, phi,dphi);
%[dcdf(:) dcdf_check(:)]
% ————————————————————————远程版本结束————————————————————————————

dcdf_check = compute_dcdf_using_phi_finite_difference(x,h,dh,g,dg, phi,dphi);
[dcdf(:) dcdf_value(:) dcdf_check(:)]


% function [predcdf,predcdf_errbnd] = compute_predcdf(h0,g0,dh0,dg0,exp_ith,...
%     PHI,DPHI,E,W,INTERVAL,pathlen,WT,EWT,NODES,...
%     DEFAULT_MAXINTERVALCOUNT,ATOL,RTOL,opstruct)
% 
% firstFunEval = true;
% 
% % subintervals
% subs = [INTERVAL(1:end-1);INTERVAL(2:end)]; 
% nsubs = size(subs,2);
% 
% % initial value
% I_sum = zeros(size(g0));
% I_error = zeros(size(g0));
% 
% while true
%     midpt = sum(subs)/2;   % midpoints of the subintervals
%     halfh = diff(subs)/2;  % half the lengths of the subintervals
%     s = NODES*halfh + midpt; % NNODES x nsubs
%     if firstFunEval
%         PHI_mat_prod = prod(PHI,1);
%         FX = imag(DPHI./PHI.*(PHI_mat_prod.*E)).*W;
%         % compute fx with kronrod weight
%         FX = reshape(FX,numel(WT),[]);
%         I_K = (WT*FX) .* halfh;
%         I_KminusG = (EWT*FX) .* halfh;
% 
%     else
%         FX = imag(DPHI./PHI.*(PHI_mat_prod.*E)).*W;
%     end
% 
% end


function [PHI,E,U,W,INTERVAL,I] = compute_pdf(phi_u,exp_ith,...
    interval,pathlen,WT,EWT,NODES,MAXINTERVALCOUNT,ATOL,RTOL,opstruct)

% subintervals
subs = [interval(1:end-1);interval(2:end)];   

% initial value 
I_sum = 0;
I_error = 0;
PHI = [];
E = [];
U = [];
W = [];
INTERVAL = interval;

while true
    midpt = sum(subs)/2;   % midpoints of the subintervals
    halfh = diff(subs)/2;  % half the lengths of the subintervals
    s = NODES*halfh + midpt; % NNODES x nsubs
    s = reshape(s,1,[]);
    tt = s ./ (1 - s);
    u = 0 + tt.^2;   % change back to original domain
    U = [U,u];
    w = 2*tt ./ (1 - s).^2;   % transform weight
    W = [W,w];

    % compute matrix phi_z
    PHI_mat = phi_u(u);
    PHI_mat_prod = prod(PHI_mat,1);

    % compute exponential part
    E_mat = exp_ith(u);

    % save PHI and E matrix
    PHI = [PHI,PHI_mat];
    E = [E,E_mat];

    % compute Real part with transform weight
    FX = real(PHI_mat_prod.*E_mat).*w;

    % compute fx with kronrod weight
    FX = reshape(FX,numel(WT),[]);
    I_K = (WT*FX) .* halfh;
    I_KminusG = (EWT*FX) .* halfh;

    % compute sum of integral
    I = sum(I_K) + I_sum;
    tol = max(ATOL,RTOL*abs(I));

    % find convergent subintervals and remove
    ndx = find(abs(I_KminusG) < (2*tol/pathlen)*abs(halfh));
    I_error = I_error + sum(I_KminusG(ndx));
    I_KminusG(ndx) = [];

    % error bound
    I_errorbnd = abs(I_error) + norm(I_KminusG,1);
    % if ~(isfinite(I) && isfinite(I_errorbnd))
    %     warning(message('NonFiniteValue'));
    %     if opstruct.ThrowOnFail
    %         error(message('integral fail'));
    %     end
    %     break
    % end
    if I_errorbnd <= tol
        break
    end

    % remove subintervals with accurate approximation
    subs(:,ndx) = [];

    % update the partial sum for the integral
    I_sum = I_sum + sum(I_K(ndx));
    midpt(ndx) = [];
    INTERVAL = unique([INTERVAL midpt]);
    if isempty(subs)
        break
    end
    % nsubs = 2*size(subs,2);
    % if nsubs > MAXINTERVALCOUNT
    %     warning(message('MaxIntervalCountReached',...
    %         sprintf('%9.1e',I_errorbnd)));
    %     if opstruct.ThrowOnFail
    %         error(message('integral fail'));
    %     end
    %     break
    % end
    subs = reshape([subs(1,:);midpt;midpt;subs(2,:)],2,[]);
end
