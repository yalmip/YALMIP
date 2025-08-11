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
        % First parameter is name, so remove that
        parameters = {distribution.parameters{2:end}};
        mixtureweights = distribution.mixture;
        
        % Evaluate terms at x
        g = funcs.g(x);
        h = funcs.h(x);
        
        % Define char. func for linear combination z = g(x)^Tw, and insert
        % the numerical parameters (means, shapes etc)  into the function
        % definition of the characteristic functions
        if ~isa(parameters{1},'cell')
            phi_z = @(t) prod(phi(g(:)*t,parameters{:}),1);        
        else
            % This is a mixture!           
             % Create mixture characteristic function...             
             phi_mixture = @(t) createMixtureSum(phi,g(:)*t,mixtureweights,parameters{:});
             phi_z = @(t) prod(phi_mixture(t),1);
        end
        
        % Compute the cdf
        varargout{1} = compute_cdf_using_phi(-h, phi_z);
        
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
        % First parameter is name, so remove that
        parameters = {distribution.parameters{2:end}};
        mixtureweights = distribution.mixture;
        % Create a function which computes gradient at x
        if isempty(mixtureweights)
            operator.derivative = @(x)compute_dcdf_using_phi(x,funcs.h,funcs.dh,funcs.g,funcs.dg,phi,dphi,parameters,mixtureweights);
        end
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};
        
    otherwise
end

function cdf = compute_cdf_using_phi(y,phi)
% Perform the inverse Fourier transform to obtain the CDF
% for Probability(z <= y) where z has characterstic function phi(t)
% Not much to exploit so for now we use built-in integrator
% Later this will be computed together with derivatives etc but this
% requires some global logic to sync those when solver wants stuff
integrand = @(t) imag(phi(t) .* exp(-1i * t * y)./t);
cdf =  .5-integral(integrand,0, inf)/pi;


function dcdf = compute_dcdf_using_phi(x,h,dh,g,dg,phi,dphi,parameters,mixtureweights)
% Compute derivative of Probability(h(x)+g(x)'*w <= 0)

% Evaluate the operators
h0  = h(x);
dh0 = dh(x);
g0  = g(x);
dg0 = dg(x);

% Create characteristic functions etc
if ~isa(parameters{1},'cell')
    phi_z   = @(t) prod(phi(g0(:)*t,parameters{:}),1);
    phi_    = @(t) phi(g0(:)*t,parameters{:});
    dphi_   = @(t) dphi(g0(:)*t,parameters{:});
    exp_ith = @(t) exp(1i*t*h0);
else
    % This is a mixture. The characteristic of the mixture has to be
    % created etc
      
end

<<<<<<< HEAD
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
=======
% Compute f_z(-h0) using Gil-Pelaez
% this will be moved to be computed together with derivative instead but
% for now we do the double work to keep code simple. Also, for now we use
% built-in integral
>>>>>>> 3722a8cddb992dccbbc8313888027791ab1dd123
integrand_pdf = @(t) real(exp_ith(t) .* phi_z(t));
if nnz(dh0)==0
    pdf_val = 0;
else
    pdf_val = (1/pi)*integral(integrand_pdf,0,inf);
end

if nnz(dg0)==0
    dcdf = (-pdf_val.*dh0');
else            
    terms = (-1/pi)*dcdfintegralEvaluator(phi_,dphi_,exp_ith);
    dcdf = (-pdf_val.*dh0') + terms'*dg0;
end

<<<<<<< HEAD
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
=======
function I = dcdfintegralEvaluator(phi_,dphi_,exp_ith)
% Compute integral on the transformed domain 0->1
% Could be various variants, but only one available now
I = adaptiveGK715(phi_,dphi_,exp_ith,0,1,setupGK715);

function I = adaptiveGK715(phi_,dphi_,exp_ith,a,b,GK715)
% We have sub-divided down to a->b. Map standard Gauss-Kronrad points to
% this interval and then map those to original domain 0->inf
center = (a+b)/2;
halfwidth = (b-a)/2;
s = GK715.Nodes*halfwidth + center; 
% and now map to original domain using nonlinear transformation
t = (s./(1-s)).^2;
coordinateChangeJacobian = (2*s./(1-s).^3).*halfwidth;

% Evaluate all functions in all points
PHI    = phi_(t);
dPHI   = dphi_(t);
PHIZ   = prod(PHI,1);
RelPHI = dPHI./PHI;
% Careful edge-case when phi is 0
RelPHI(isnan(RelPHI))=0;

% Product characteristic, exponential and transformation can be combined
Weight = (PHIZ.*exp_ith(t).*coordinateChangeJacobian);
% Integral = sum of weighted columns
I_15 = imag(RelPHI*(Weight.*GK715.Weights_15).');
I_7 = imag(RelPHI(:,2:2:end)*(Weight(2:2:end).*GK715.Weights_7).');
error = abs(I_15-I_7);
% Simple stopping criteria for now, and we sub-divide in all despite some
% might be done already or regions being all 0. Will be optimized later
if all(error <= (b-a)*1e-6)
    I = I_15;    
else
    I1 = adaptiveGK715(phi_,dphi_,exp_ith,a,(a+b)/2,GK715);
    I2 = adaptiveGK715(phi_,dphi_,exp_ith,(a+b)/2,b,GK715);
    I = I1+I2;    
end

function GK715 = setupGK715
GK715.Nodes = [ ...
    -0.9914553711208126, -0.9491079123427585, -0.8648644233597691, ...
    -0.7415311855993944, -0.5860872354676911, -0.4058451513773972, ...
    -0.2077849550078985, 0, 0.2077849550078985, ...
    0.4058451513773972, 0.5860872354676911, 0.7415311855993944, ...
    0.8648644233597691, 0.9491079123427585, 0.9914553711208126];
GK715.Weights_15 = [ ...
    0.02293532201052922, 0.06309209262997855, 0.1047900103222502, ...
    0.1406532597155259, 0.1690047266392679, 0.1903505780647854, ...
    0.2044329400752989, 0.2094821410847278, 0.2044329400752989, ...
    0.1903505780647854, 0.1690047266392679, 0.1406532597155259, ...
    0.1047900103222502, 0.06309209262997855, 0.02293532201052922];
GK715.Weights_7 = [ ...
    0.1294849661688697, 0.2797053914892767, 0.3818300505051189, ...
    0.4179591836734694, 0.3818300505051189, 0.2797053914892767, ...
    0.1294849661688697];


function C = createMixtureSum(f, g_t,mixtureweights,varargin)
% Characteristic functions sum_j alpha_j phi(g*t,parameters_j)
C = [];
% varargin are all distribution parameters, each containing a cell with
% the indvidual component parameters
% Every element can have different mixture weights, but the number of
% components is the same for all elements
n_mixtures = length(varargin{1});
C = 0;
for j = 1:n_mixtures
    % Extract parameters for this mixture component
    jthParameters = cellfun(@(c) c{j}, varargin, 'UniformOutput', false);
    % Evaluate characteristic function. Note that g(x)*t can be a row
    % vector as an integral evaluator might evaluate multiple points
    C = C + mixtureweights(:,j).*f(g_t,jthParameters{:});
end
>>>>>>> 3722a8cddb992dccbbc8313888027791ab1dd123
