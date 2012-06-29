function  data = datasaver(varargin)

persistent saved

if nargin > 3
    saved.savedfunctiondata = varargin{1};
    saved.savedfunctionGradientdata = varargin{2};
    saved.savedfunctionHessiandata = varargin{3};
    saved.savedconstraintdata = varargin{4};
    saved.savedconstraintGradientdata = varargin{5};
    saved.savedconstraintHessiandata = varargin{6};
    saved.m = size(saved.savedconstraintdata.F_struc,1);
    saved.n = length(saved.savedfunctiondata.linearindicies);
    saved.monomials = [];
    saved.x = randn(saved.n,1);
    return
end

x = varargin{2};
if ~isequal(x,saved.x)
    % Recompute the monomial terms for new x
    saved.monomials=ones(1,size(saved.savedfunctiondata.monomtable,2));
    saved.monomials(saved.savedfunctiondata.linearindicies) = x(:)';
   % changed = find(x-saved.x);
   % changed = saved.savedfunctiondata.linearindicies(changed);
   % changedmonomials = find(any(saved.savedfunctiondata.monomtable(saved.savedfunctiondata.nonlinearindicies,changed),2));
    
    for i = 1:length(saved.savedfunctiondata.nonlinearindicies)
        saved.monomials(saved.savedfunctiondata.nonlinearindicies(i)) = prod(saved.monomials.^saved.savedfunctiondata.monomtable(saved.savedfunctiondata.nonlinearindicies(i),:));
    end
    saved.x = x;
end

switch varargin{1}
    case 1
        F = saved.savedfunctiondata.F_struc;
        x = varargin{2};
    case 2
        F = saved.savedfunctionGradientdata.F_struc;
    case 3
        F = saved.savedfunctionHessiandata.F_struc;
    case 4
        if nargin == 2
            % ipopt call
            F = saved.savedconstraintdata.F_struc;
        else
            % pennlp call
            F = saved.savedconstraintdata.F_struc(1+varargin{3},:);
        end
    case 5
        F = saved.savedconstraintGradientdata.F_struc;
        start = 1+(varargin{3}-1)*saved.n;
        F = F(varargin{3}:saved.m:end,:);
    case 6
        F = saved.savedconstraintHessiandata.F_struc;
        F = full(F((varargin{3}-1)+(1:saved.m:saved.m*saved.n^2),:));        
    otherwise
end

data = (F*[1;saved.monomials(:)]);