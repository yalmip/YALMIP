function Intervals = redistribute_grid(Intervals, Errors, alpha)
% Redistribute grid points based on error equidistribution
% Intervals: 2 x N array of interval endpoints
% Errors: M x N array of error vectors
% alpha: damping factor (0 < alpha <= 1)

if nargin < 3
    alpha = 0.1;
end

if isempty(Intervals) || size(Intervals,2) < 2
    return;
end

% Extract current grid points
current_grid = [Intervals(1, :) Intervals(2, end)];

% Calculate scalar error metric for each interval
if exist('Errors', 'var') && ~isempty(Errors)
    m = mean(Errors, 1);
else
    m = ones(1, size(Intervals, 2));
end

% Avoid division by zero and extreme scaling
m = max(m, 1e-12);

% We want m_i * h_i = C
inv_m = 1 ./ m;
C = 1 / sum(inv_m);
fair_lengths = C * inv_m;

% Construct fair grid
fair_grid = [0 cumsum(fair_lengths)];
fair_grid(end) = 1; % Ensure exact 1

% Damped update
new_grid = current_grid + alpha * (fair_grid - current_grid);

% Reconstruct intervals
Intervals = [new_grid(1:end-1); new_grid(2:end)];
end
