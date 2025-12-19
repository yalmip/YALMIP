function [MergedIntervals, MergedErrors] = merge_adaptive_intervals(Intervals, Errors, mergeStrategy)

% Errors is assumed to be M x N (M error terms, N intervals)
% Intervals is 2 x N

% Define operators for error handling
smashOperator = @(e1, e2) 10*(e1 + e2);
scalarizeOperator = @(e) max(abs(e),[], 1);

if mergeStrategy == 1
    [MergedIntervals, MergedErrors] = mergeIntervalsLeftToRight(Intervals, Errors, smashOperator, scalarizeOperator);
elseif mergeStrategy == 2
    [MergedIntervals, MergedErrors] = mergeIntervalsRightToLeft(Intervals, Errors, smashOperator, scalarizeOperator);
elseif mergeStrategy == 4
    [MergedIntervals, MergedErrors] = mergeIntervalsTreeConservative(Intervals, Errors);
else
    [MergedIntervals, MergedErrors] = mergeIntervalsMinError(Intervals, Errors, smashOperator, scalarizeOperator);
end

function [D_merged, E_merged] = mergeIntervalsLeftToRight(D, E, smashOp, scalarOp)
if size(D,2) <= 4
    D_merged = D;
    E_merged = E;
    return;
end

D_merged = D(:,1);
E_merged = E(:,1);

for i = 2:size(D,2)
    last_interval = D_merged(:,end);
    last_error = E_merged(:,end);
    
    current_interval = D(:,i);
    current_error = E(:,i);
    
    candidate_error = smashOp(last_error, current_error);
    candidate_len = current_interval(2) - last_interval(1);
    tol = max(1e-6, candidate_len*1e-6);
    
    % Check if merging would violate the minimum interval count constraint
    remaining_intervals = size(D, 2) - i;
    current_count = size(D_merged, 2);
    can_merge = (current_count + remaining_intervals) >= 4;
    
    if can_merge && scalarOp(candidate_error) <= tol
        % Merge
        D_merged(2,end) = current_interval(2);
        E_merged(:,end) = candidate_error;
    else
        % Append
        D_merged = [D_merged current_interval];
        E_merged = [E_merged current_error];
    end
end

if size(D_merged,2) < size(D,2)
    [D_merged, E_merged] = mergeIntervalsLeftToRight(D_merged, E_merged, smashOp, scalarOp);
end

function [D_merged, E_merged] = mergeIntervalsRightToLeft(D, E, smashOp, scalarOp)
if size(D,2) <= 3
    D_merged = D;
    E_merged = E;
    return;
end

D_merged = D(:,end);
E_merged = E(:,end);

for i = size(D,2)-1:-1:1
    first_merged_interval = D_merged(:,1);
    first_merged_error = E_merged(:,1);
    
    current_interval = D(:,i);
    current_error = E(:,i);
    
    candidate_error = smashOp(first_merged_error, current_error);
    candidate_len = first_merged_interval(2) - current_interval(1);
    tol = max(1e-6, candidate_len*1e-6)/2;
    
    % Check if merging would violate the minimum interval count constraint
    remaining_intervals = i - 1;
    current_count = size(D_merged, 2);
    can_merge = (current_count + remaining_intervals) >= 3;
    
    if can_merge && scalarOp(candidate_error) <= tol
        % Merge
        D_merged(1,1) = current_interval(1);
        E_merged(:,1) = candidate_error;
    else
        % Prepend
        D_merged = [current_interval D_merged];
        E_merged = [current_error E_merged];
    end
end

if size(D_merged,2) < size(D,2)
    [D_merged, E_merged] = mergeIntervalsRightToLeft(D_merged, E_merged, smashOp, scalarOp);
end

function [D, E] = mergeIntervalsMinError(D, E, smashOp, scalarOp)
if isempty(D)
    return;
end

while size(D,2) > 3
    % Compute candidate errors for all adjacent pairs
    % E is M x N, so we need to sum adjacent columns
    cand_E = smashOp(E(:,1:end-1), E(:,2:end));
    
    % Compute candidate lengths
    cand_L = D(2, 2:end) - D(1, 1:end-1);
    
    % Compute tolerances
    tols = max(1e-6, cand_L * 1e-6);
    
    % Scalarize errors for comparison (returns row vector)
    scalar_cand_E = scalarOp(cand_E);
    
    % Find valid merges
    valid_mask = (scalar_cand_E <= tols)/10;
    
    if ~any(valid_mask)
        break;
    end
    
    % Find the best merge (minimum error density) among valid ones
    valid_indices = find(valid_mask);
    scores = scalar_cand_E(valid_indices) ./ cand_L(valid_indices);
    [~, best_idx_in_valid] = min(scores);
    best_idx = valid_indices(best_idx_in_valid);
    
    % Perform merge
    D(2, best_idx) = D(2, best_idx+1);
    E(:, best_idx) = cand_E(:, best_idx);
    
    % Remove the merged interval
    D(:, best_idx+1) = [];
    E(:, best_idx+1) = [];
end

function [D, E] = mergeIntervalsTreeConservative(D, E)
% Merges siblings (adjacent intervals of same length) if their combined error
% is low enough. This preserves the tree structure of the adaptive grid.

if isempty(D)
    return;
end

% Safety factor C < 1
C = 0.5;

didMerge = true;
while didMerge && size(D,2) > 3
    didMerge = false;
    i = 1;
    while i < size(D,2)
        % Check if siblings (approx same length)
        len1 = D(2,i) - D(1,i);
        len2 = D(2,i+1) - D(1,i+1);
        
        if 1%abs(len1 - len2) < 1e-1
            % Candidate for merge
            cand_E = E(:,i) + E(:,i+1);
            cand_L = len1 + len2;
            
            % Tolerance for this merged interval
            tol = max(1e-6, cand_L * 1e-6);
            
            % Vector check: All components must be safe
            % (E_L + E_R) < C * tol
            if all(cand_E < .1 * tol)
                % Merge
                D(2,i) = D(2,i+1);
                E(:,i) = cand_E;
                
                % Remove i+1
                D(:,i+1) = [];
                E(:,i+1) = [];
                
                didMerge = true;
                
                % Check constraint
                if size(D,2) <= 3
                    return;
                end
            else
                i = i + 1;
            end
        else
            i = i + 1;
        end
    end
end
