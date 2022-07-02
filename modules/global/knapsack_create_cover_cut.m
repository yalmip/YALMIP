function cut = knapsack_create_cover_cut(a,b,x,alg,gubs)
% Derive cover inequalities cut*[1;x]>=0 for a*x <= b

if (all(a>=0) && sum(a)<=b) || (all(a<=0) && (b>=0))
    % quick exit sum x <= large or sum x >= neg
    cut = [];
    return
end

if nargin < 3
    x = ones(length(a),1);
    alg = 'crowder';
    gubs = [];
elseif nargin < 4
    alg = 'crowder';
    gubs = [];
elseif nargin < 5
    gubs = spalloc(1,length(a),0);
end

if ~all(a)
    nz = find(a);
    a_ = a(nz); 
    cut_ = knapsack_create_cover_cut(a_,b,x(nz),alg,gubs(nz));
    if ~isempty(cut_)
        cut = spalloc(1,length(a)+1,0);
        cut(1) = cut_(1);
        cut(1+nz) = cut_(2:end);
    else
        cut = [];
    end
    return
elseif any(a<0)
    % Negative values in a*x <= b
    % models -a*x >= -b i.e. typically something like sum stuff >= bound
    %        -a*(1-y) >= -b  
    %        -a*y <= b-sum(a)
    neg = find(a < 0);
    a_ = a;
    a_(neg) = -a(neg);
    x_ = x;
    x_(neg) = 1-x(neg);
    cut = knapsack_create_cover_cut(a_,b-sum(a(neg)),x_,alg,gubs);
    % we know have a cut, cut(1) + cut(2:end)*y >=0
    %                     cut(1) + cut(2:end)(1-x) >= 0
    %                     cut(1) + sum(cut(2:end)) + (-cut(2))*x >= 0
    if ~isempty(cut)
        cut(1) = cut(1)+sum(cut(1 + neg));
        cut(1 + neg) = -cut(1 + neg);
    end
    return
end
    
switch alg
    case {'crowder',''}
        % Heuristics from H. Crowder, E. Johnson, M. Padberg
        % Solving large-scale 0–1 linear programming programs
        [val,loc] = sort((1-x)./a(:),'ascend');
    case 'gu'
        [val,loc] = sort((x),'descend');
    case 'gu-reverse'
        [val,loc] = sort((x),'ascend');
                
    otherwise
        error('Unsupport cover cut separation')
end
% Initial cover
C = min(find(cumsum(a(loc))>b));

% Apply Balas lifting
q = knapsack_cover_lift_balas(a,loc(1:C));
% qL = knapsack_cover_lift_letchford(a,b,loc(1:C));
% if q~=qL
%    'HEJ'
% end
% Return row where row*[1;x] hopefully is violated
cut = spalloc(1,length(a)+1,0);
cut(1) = C-1;
cut(2:end)=-q;




