function [g,dg] = sdpLayerProjectorOracle(X,dX,algorithm,nConstraints,clusterTolerance,signedDistancePowers)
% Values and derivatives for invariant spectral SDP scalarizations.
%
% dX contains vectorized matrix derivatives, one derivative per column.
% dg is returned with one constraint per row.

X = full((X + X')/2);
[V,D] = eig(X);
[eigenvalues,order] = sort(real(diag(D)),'ascend');
V = real(V(:,order));

scale = max(1,max(abs(eigenvalues)));
clusters = determineClusters(eigenvalues,clusterTolerance*scale);

if nargin < 6
    signedDistancePowers = [1 2 inf];
end

switch lower(algorithm)
    case 'projector'
        active = clusters{1};
        g = -eigenvalues(1);
        if nargout > 1
            dg = -projectorGradient(V(:,active),dX);
        end

    case 'signed-distance'
        negative = find(eigenvalues < 0);
        if isempty(negative)
            active = clusters{1};
            g = -eigenvalues(1);
            if nargout > 1
                dg = -projectorGradient(V(:,active),dX);
            end
        else
            negativeEigenvalues = eigenvalues(negative);
            g = norm(negativeEigenvalues);
            if nargout > 1
                weights = negativeEigenvalues/g;
                dg = weightedSpectralGradient(V(:,negative),weights,dX);
            end
        end

    case 'signed-distance-multi'
        powers = validateSignedDistancePowers(signedDistancePowers);
        negative = find(eigenvalues < 0);
        g = zeros(length(powers),1);
        if nargout > 1
            dg = zeros(length(powers),size(dX,2));
        end
        if isempty(negative)
            active = clusters{1};
            g(:) = -eigenvalues(1);
            if nargout > 1
                row = -projectorGradient(V(:,active),dX);
                dg(:,:) = repmat(row,length(powers),1);
            end
        else
            violations = -eigenvalues(negative);
            for i = 1:length(powers)
                p = powers(i);
                if isinf(p)
                    active = clusters{1};
                    g(i) = -eigenvalues(1);
                    if nargout > 1
                        dg(i,:) = -projectorGradient(V(:,active),dX);
                    end
                else
                    g(i) = norm(violations,p);
                    if nargout > 1
                        weights = -(violations/g(i)).^(p-1);
                        dg(i,:) = weightedSpectralGradient(V(:,negative),weights,dX);
                    end
                end
            end
        end

    case 'projector-multi'
        nConstraints = min(nConstraints,length(eigenvalues));
        g = zeros(nConstraints,1);
        if nargout > 1
            dg = zeros(nConstraints,size(dX,2));
        end
        for i = 1:length(clusters)
            members = clusters{i};
            slots = members(members <= nConstraints);
            if ~isempty(slots)
                g(slots) = -mean(eigenvalues(members));
                if nargout > 1
                    row = -projectorGradient(V(:,members),dX);
                    dg(slots,:) = repmat(row,length(slots),1);
                end
            end
        end

    otherwise
        error('Unknown invariant spectral SDP layer algorithm ''%s''.',algorithm);
end

function row = projectorGradient(V,dX)
n = size(V,1);
row = zeros(1,size(dX,2));
for j = 1:size(dX,2)
    Xj = reshape(dX(:,j),n,n);
    row(j) = trace(V'*Xj*V)/size(V,2);
end

function row = weightedSpectralGradient(V,weights,dX)
n = size(V,1);
row = zeros(1,size(dX,2));
for j = 1:size(dX,2)
    Xj = reshape(dX(:,j),n,n);
    row(j) = weights'*diag(V'*Xj*V);
end

function powers = validateSignedDistancePowers(powers)
powers = powers(:)';
if isempty(powers) || ~isreal(powers) || any(isnan(powers)) || any(powers < 1)
    error('slayer.p must contain real exponents in [1,inf].');
end

function clusters = determineClusters(eigenvalues,tolerance)
clusters = {1};
current = 1;
for i = 2:length(eigenvalues)
    if abs(eigenvalues(i)-eigenvalues(i-1)) <= tolerance
        clusters{current}(end+1) = i;
    else
        current = current + 1;
        clusters{current} = i;
    end
end
