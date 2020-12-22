function [W, y] = csp(S1, S2, varargin)
% Compute spatial filters from covariance matrices
%
% [W, y] = csp(S1, S2, 'name', value)
%   It computes the spatial filters W from covariance matrices S1 and S2.
%   S1 and S2 can be also absolute paths to the covariance matrices. It
%   uses 'name', value parameters to specify the CSP method and other
%   optional parameters like the number of spatial filters.
%
% Parameters
% ----------
% S1 (double, str):
%   Average covariance matrix of condition 1, with size [C C], where C is 
%   the number of channels. It can be also the absolute path to a matrix S
%   that has size [C C N], with N being the number of trials or windows.
% S2 (double, str):
%   Average covariance matrix of condition 2, with size [C C], where C is 
%   the number of channels. It can be also the absolute path to a matrix S
%   that has size [C C N], with N being the number of trials or windows.
%
% Optional name-value parameters
% ------------------------------
% 'method': 'regular', 'max-sliced-Bures' (str)
%   Name of method used to compute W. Currently supported:
%       'regular': Maximization of E1k/E2k, with Eck being the variance
%       (energy) of the projected (filtered) signal under condition c, with
%       c = {1,2}. Columns of W are sorted in descending order of d1, the 
%       eigenvalues of S1 (the average covariance matrix of condition 1).
%       d1 + d2 = 1.
%       'max-sliced-Bures': Finds w1 and w2 that maximizes g(w1)-h(w1), and
%       h(w2)-g(w2), respectively, with g(w) = sqrt(w'*S1*w) and 
%       h(w) = sqrt(w'*S2*w). k must be equal to 2.
%       'betadiv': Uses the symmetric beta divergence to estimate the
%       spatial filters. See "Robust Spatial Filtering with Beta 
%       Divergence", by Samek et al., 2013.
% 'k': int
%   Number of spatial filters requested. If the method is 'regular' or 
%   'max-sliced-Bures', k must be even. Default: return all.
% 'beta': float
%   Beta parameter in the beta-divCSP method.
% 'wfname': [], or wfname (str)
%   If [], filtcov will not save W to disk. If wfname is provided, it will
%   save it as wfname.
%
% Returns
% ------------------------
% W (double):
%   A matrix of size [C k], with C the number of channels and k the
%   number of spatial filters. If k is not specified, k = C.
%

%% Defaults for optional parameters
method = 'regular'; k = []; beta = 0.5; wfname = [];

%% Parse arguments
if nargin < 2
    error('Wrong number of arguments');
elseif nargin > 2
    n_varargin = length(varargin);
    assert(mod(n_varargin,2) == 0, 'Wrong number of optional arguments');
    for i_opt = 1:2:n_varargin
        switch varargin{i_opt}
            case 'method'
                method = varargin{i_opt+1};
            case 'k'
                k = varargin{i_opt+1};                
            case 'beta'
                beta = varargin{i_opt+1};
            case 'wfname'
                wfname = varargin{i_opt+1};
            otherwise
                error('Wrong optional arguments');
        end
    end
end

switch method
    case 'regular'
        [S1, S2] = checkS(S1, S2);
        assert(mod(k,2)==0, 'k must be even');
        [W, y] = vanilla_csp(S1, S2, 'k', k);
    case 'max-sliced-Bures'
        assert(k==2, 'Max-sliced Bures only returns two CSP filters');
        [S1, S2] = checkS(S1, S2);
        [W, y] = max_sliced_Bures(S1, S2);
    case 'betadiv'
        [W, ~, y] = betadivcsp(S1, S2, beta, 'k', k, 'quiet', false);
end

if ~isempty(wfname)
    save(wfname, 'W');
end

end

function [S1, S2] = checkS(S1, S2)
if ischar(S1) && ischar(S2)    
        Saux = load(S1, 'S');
        S1 = mean(Saux.S,3);
        Saux = load(S2, 'S');
        S2 = mean(Saux.S,3);
end
end