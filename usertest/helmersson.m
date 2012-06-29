function pass = helmersson

% this example is from scherer, lpv control and full block
% multipliers, automatica 37 (2001), pp. 361-375.

clear all

if ~exist ('a', 'var'), 
  a = 1.0; 
elseif size (a) ~= [1 1], 
  a = 1.0; 
end

sz = 3;
if ~exist ('ndd', 'var'), 
  ndd = sz; 
elseif size (ndd) ~= [1 1], 
  ndd = sz; 
end

if ~exist ('nddr', 'var'), 
  nddr = sz; 
elseif size (nddr) ~= [1 1], 
  nddr = sz; 
end

ddr = 2.0/nddr;
dd = ddr/ndd; 

if ~exist('slack_factor', 'var'), 
  slack_factor = 1.0; 
elseif size (slack_factor) ~= [1 1], 
  slack_factor = 1.0; 
end

A = -1;
B = [1 1 1 1 0 1];
C = [0; 0; 0; 0; 0; 1];
D = [0     1 0   1 0 0;
     0.5   0 0.5 0 1 0;
     2*a   0 a   0 0 0;
     0  -2*a 0  -a 0 1;
     0     1  0  0 0 0;
     1     0  0  0 0 0];

I = eye(size([A B],2));
     
blk = [1 4 1 1];  % states deltas perf control
ABCD = [A B; C D];


deltar = 0.9*(-1+0.5*ddr:ddr:1-0.5*ddr);

Deltars = zeros (4, 4, length(deltar)^2);
m = 0;
nd = length (deltar);
for i1 = 1:nd,
  d1 = deltar(i1);
  for i2 = 1:nd,
    d2 = deltar(i2);
    m = m + 1;
    Deltars(:,:,m) = blkdiag (d1*eye(2), d2*eye(2));
  end
end

nr = size(Deltars,3);

delta = 0.9*(-0.5*ddr:dd:0.5*ddr);

maxsize = max ([0 diff(delta); diff(delta) 0]);
Deltas = zeros (4, 4, length(delta)^2, nr);
slack = slack_factor * ones (length(delta)^2, 2);

for ir = 1:nr
  maxsize = max ([0 diff(delta); diff(delta) 0]);
  slack = slack_factor * ones (length(delta)^2, 2);

  m = 1;
  nd = length (delta);
  for i1 = 1:nd,
    d1 = delta(i1);
    for i2 = 1:nd,
      d2 = delta(i2);
      Deltas(:,:,m,ir) = blkdiag (d1*eye(2), d2*eye(2)) + Deltars(:,:,ir);
      slack(m,:) = slack_factor * [maxsize(i1) maxsize(i2)];
      m = m + 1;
    end
  end

  DDeltas = [];
  DDeltas(:,:,1) = blkdiag (eye(2), 0*eye(2));
  DDeltas(:,:,2) = blkdiag (0*eye(2), eye(2));
end
  
[P, R, gam, F] = synlpv (ABCD, blk, Deltas, DDeltas, slack);

% evaluation

Pd = P(3:10,3:10);
Rd = R(3:10,3:10);

pass = norm(gam-1.527)<1e-2;


function [P_feas, R_feas, gam_feas, F] = synlpv (ABCD, blk, Deltas, ...
					      DDeltas, slack)

if nargin < 3,
  disp ('usage:  [P, R, gam] = synlpv (A, blk, Deltas, DDeltas, slack)');
  return;
end

if nargin < 4, DDeltas = []; end

n = sum(blk(1:3));
A = ABCD(1:n,1:n);
B = ABCD(1:n,n+1:end);
C = ABCD(n+1:end,1:n);

n1 = blk(1);  % states
n2 = blk(2);  % deltas
n3 = blk(3);  % H-inf
[n, k]      = size (A);
[kd, nd, q, r] = size (Deltas);
[kdd, ndd, qd] = size (DDeltas);
if isempty (DDeltas), qd = 0; end

if n ~= n1+n2+n3, error ('blk not compatible with A'); end
if n ~= n1+n2+n3, error ('blk not compatible with A'); end
if k ~= n1+n2+n3, error ('blk not compatible with A'); end

if n2 ~= nd | n2 ~= kd, error ('blk not compatible with Deltas'); end

if qd > 0 & (kdd ~= kd | ndd ~= nd), 
  error ('DDeltas not compatible with Deltas');
end

if nargin < 5 & qd > 1,
  slack = ones(q,qd);
end

if length (slack) ~= q, 
  error ('slack not compatible with Deltas');
end

slackvar = any (slack(:) > 0);

ns = n1;
ks = n1;   % state block is always square

np = n3;   % performance block (H-inf)
kp = n3;

tic
yalmip ('clear');

% define the block diagonal structure of P and R

Ps = sdpvar (ns, ns, 'symmetric');
Rs = sdpvar (ns, ns, 'symmetric');
zs = zeros(ns,ns);

gam  = sdpvar (1, 1);
gami = sdpvar (1, 1);

for ir = 1:r, 
  Pd(ir).var = sdpvar (nd+kd, nd+kd, 'symmetric');
  P(ir).var  = blkdiag ([zs Ps; Ps zs], Pd(ir).var);
  Rd(ir).var = sdpvar (nd+kd, nd+kd, 'symmetric');
  R(ir).var  = blkdiag ([zs Rs; Rs zs], Rd(ir).var);

  % Hinf outputs and inputs
  P(ir).var = blkdiag (P(ir).var, gami*eye(np), -gam*eye(kp));
  R(ir).var = blkdiag (R(ir).var, gam*eye(np), -gami*eye(kp)); 

  if slackvar,
    % define slack variables of Delta inequalities
    for iq = 1:qd,
      DDP(iq,ir).var = sdpvar (nd, nd, 'symmetric');
      DDR(iq,ir).var = sdpvar (nd, nd, 'symmetric');
    end
  end
end

Pswap = [1:ns n+(1:ks) ns+(1:nd) n+ns+(1:kd) ...
	ns+nd+(1:np) n+ns+nd+(1:kp)];

Rswap = [1:ks k+(1:ns) ks+(1:kd) k+ns+(1:nd) ...
	ks+kd+(1:kp) k+ns+nd+(1:np)];

AI = [A; eye(k)];
IA = [eye(n); -A'];

ACp = AI(Pswap,:)*null(C);

ABp = IA(Rswap,:)*null(B');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gam*gami > 1
F = lmi;

% [gam 1; 1 gami] > 0, 
% which is equivalent to gam > 0 and gami > 1/gam
% even if gami > 1/gam we can replace it by 1/gam and still
% satisfying the quadratic constraints since gami terms appears as
% in the P and R inequalities as
% ... + U' gami U < 0 and ... - V' gami V > 0, respectively
F = addlmi (F, [gam 1; 1 gami]);

% state-space multipliers defined by [Ps I; I Rs] > 0
F = addlmi (F, [Ps eye(ns); eye(ns) Rs]);

for ir = 1:r
  % Cp' [A;I]' P(ir) [A;I] Cp < 0
  F = addlmi (F, -ACp'*P(ir).var*ACp);

  % Bp' [I;A']' R(ir) [I;A'] Cp > 0
  F = addlmi (F, ABp'*R(ir).var*ABp);

  if slackvar,
    % slack variables
    for iq = 1:qd,
      % DDP(iq,ir) > [0; DD]' Pd(ir) [0; DD]
      ZDD = [zeros(nd,nd); DDeltas(:,:,iq)];

      F = addlmi (F, DDP(iq,ir).var - ZDD'*Pd(ir).var*ZDD);

      % DDP(iq,ir) > 0
      F = addlmi (F, DDP(iq,ir).var);
  
      % DDR(iq,ir) + [-DD'; 0]' Rd(ir) [-DD'; 0] > 0
      DDZ = [-DDeltas(:,:,iq)'; zeros(kd,kd)];
      F = addlmi (F, DDR(iq,ir).var + DDZ'*Rd(ir).var*DDZ);

      % DDR(iq,ir) > 0
      F = addlmi (F, DDR(iq,ir).var);
    end
  end

  for iq = 1:q,
    % [I; D]' Pd(ir) [I; D] > 0.25 * sum DDP
    sumDDP = 0;
    if slackvar,
      for iqd = 1:qd,
        sumDDP = sumDDP + slack(iq,iqd)^2*DDP(iqd,ir).var;
      end
    end
    
    ID = [eye(nd); Deltas(:,:,iq,ir)];
    F = addlmi (F, ID'*Pd(ir).var*ID - 0.25*sumDDP);

    % [-D'; I]' Rd(ir) [-D'; I] + 0.25 * sum DDR < 0
    sumDDR = 0;
    if slackvar,
      for iqd = 1:qd,
        sumDDR = sumDDR + slack(iq,iqd)^2*DDR(iqd,ir).var;
      end
    end
    DI = [-Deltas(:,:,iq,ir)'; eye(kd)];
    F = addlmi (F, -DI'*Rd(ir).var*DI - 0.25*sumDDR);
  end
end

toc

% minimize gam
%[c,A,b,K] = sedumidata(F,gam)
%save dummy c A b K
%break
sol = solvesdp (F,gam)

if isempty (sol),
  disp ('no solution found');
  P = [];
  R = [];
else
  gam_feas = sdp2mat (gam, sol);  
  gami_feas = sdp2mat (gami, sol);
  disp (sprintf ('gam = %8.3e', gam_feas))
  disp (sprintf ('gam*gami - 1 = %8.3e', gam_feas*gami_feas - 1))
  P_feas = sdp2mat (P(1).var, sol);
  R_feas = sdp2mat (R(1).var, sol);
end



