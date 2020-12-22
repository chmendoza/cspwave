function [W, y] = vanilla_csp(S1, S2, varargin)
% See "Optimizing Spatial filters for Robust EEG Single-Trial
% Analysis"

k = [];
if nargin == 4 % 'k',k
    k = varargin{2};
end

fprintf('Computing CSP filters...\n');

n_chan = size(S1, 1);

%% Solve general eigenvalue problem
St = S1 + S2;
[W, D] = eig(S1, St);
[~, isort] = sort(diag(D), 'descend');
W = W(:,isort);

D1 = W' * S1 * W; D2 = W' * S2 * W;
assert(isdiag(D1)); assert(isdiag(D2));
D1 = diag(D1);  D2 = diag(D2);
isClose2One = abs(D1+D2-ones(n_chan,1)) < 1e-10;
assert(all(isClose2One));
gtind = D1 >= D2;
y = [D1(gtind)./D2(gtind); D2(~gtind)./D1(~gtind)];

if ~isempty(k)
    k = k/2;
    W = [W(:,1:k) W(:,end-k+1:end)];
end
end

function tf = isdiag(X)
% Relaxed version of built-in isdiag()
X = X - diag(diag(X));
if max(abs(X(:))) > 1e-10
    tf = false;
else
    tf = true;
end
end