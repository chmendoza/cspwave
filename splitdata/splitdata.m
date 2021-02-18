function [fnames, i_start] = splitdata(dirpath, winlen, n_train, n_test, varargin)
% Pick multichannel windows of length winlen at random for training and testing
% of CSP method
%
% [fnames, i_start] = splitdata(dirpath, winlen, n_train, n_test, varargin)
%
% Parameters
% ----------
% dirpath (str):
%   Absolute path to the data folder
% winlen (int):
%   Length of the instance (window)
% n_train (int):
%   Number of training instances
% n_test (int):
%   Number of test instances
% 
% Optional parameters
% -------------------
% overlap (logical):
%   true if overlapping windows for train/test set are allowed, false
%   otherwise. Default: false.
% n_subdiv (int):
%   Number of non-overlapping subintervals where to pick from the
%   overlapping train/test windows. The length of the subintervals on each 
%   epoch is equal. Subintervals are assigned at random to train/test sets.
%   It must be an even number. Default: 2.
% start_gap (int):
%   Gap left at the beginning of each epoch to pick a random start. Default
%   is equal to winlen.
%
% Returns
% -------
% fnames (str): cell array with the filenames of the epochs from where 
%               windows were picked
%
% i_start (uint32): cell array with start index of the windows within each epoch
%      
% i_start{i,k}(j) is the start index of the j-th window in the file
% fnames{i}, with k=1 for train set and k=2 for test set
%
% The folder in dirpath contains files with the pattern name
% rx<epoch_id>.mat, where epoch_id is an integer. Each file has the
% variable 'epoch', an CxT matrix, with C the number of channels, and T
% the number of time points. T might vary among epochs (files).
%
% If overlap is false, the number of random starts (e.g., Monte Carlo 
% simulations) is limited to winlen-1. There are much more random starts if
% overlap is true.

%% Defaults for optional parameters
overlap = false; n_subdiv = 2; start_gap = winlen;

%% Parse arguments
if nargin < 4
    error('Wrong number of arguments');
elseif nargin > 4
    n_varargin = length(varargin);
    assert(mod(n_varargin,2) == 0, 'Wrong number of optional arguments');
    for i_opt = 1:2:n_varargin
        switch varargin{i_opt}
            case 'overlap'
                overlap = varargin{i_opt+1};
            case 'n_subdiv'
                n_subdiv = varargin{i_opt+1};
            case 'start_gap'
                start_gap = varargin{i_opt+1};            
            otherwise
                error('Wrong optional arguments');
        end
    end
end
         
%% Get file list
expr = '^rx(?<epoch_id>\d+).mat$';
S = dir(dirpath);
fnames = {S.name};
epoch_files = ~cellfun(@isempty, regexp(fnames, expr));
fnames = fnames(epoch_files);
n_epochs = length(fnames);

%% Pick number of epochs (files)
% If the number of epochs is too big, pick only max(n_train, n_test). In
% that way, at least one window is picked per epoch
max_epochs = max(n_test, n_train);
if n_epochs > max_epochs
    % Randomly shuffle the file names
    fnames = fnames(randperm(n_epochs, max_epochs));
    n_epochs = max_epochs;
end

%% Get the length of each epoch to estimate how many windows should be picked
pnts = zeros(n_epochs,1);
for i_epoch = 1:n_epochs
    %% Load the data    
    fpath = fullfile(dirpath, fnames{i_epoch});
    mObj = matfile(fpath);
    pnts(i_epoch) = size(mObj, 'epoch', 2);    
end
[pnts, isort] = sort(pnts, 'descend');
fnames = fnames(isort);

% Discard epochs that can not accomodate at least 1 window/subdiv or 1
% window
if overlap
    % FIXME: This might be wrong
    ind = pnts < (winlen+start_gap)*n_subdiv;
else
    ind = pnts < winlen+start_gap;
end
fprintf('Discarding %d out %d epochs that are not long enoug\n', ...
    sum(ind), n_epochs);
pnts(ind) = [];
fnames(ind) = [];

%% Make sure that there is enough data available
L = sum(pnts);  % Total number of samples available (per channel)
if overlap
    n_win = n_train + n_test;
    % Start indexes can go in [1 l-start_gap], where l is the length of a
    % subdivision of an epoch. So, start_gap samples are reserved per
    % subdivision per epoch.
    % FIXME: how big can be winlen?
    max_win = L - n_subdiv * start_gap * n_epochs;
    data_prctage = n_win / max_win;
else
    % add 1 start_gap/epoch for a random start. NOTE: This limits the number of
    % random starts (Monte Carlos simulations) to start_gap - 1:
    n_win = n_train + n_test;
    data_prctage = winlen * n_win / (L-start_gap*n_epochs);  % percentage of data requested
end

assert(data_prctage<=1, 'Not enough data in %d epochs!', n_epochs);

if overlap
    tot_n_win = floor((pnts-n_subdiv*start_gap)*data_prctage);
else
    tot_n_win = floor((pnts-start_gap)*data_prctage/winlen);
end

train_prctage = n_train / (n_train + n_test);

n_win = zeros(n_epochs, 2);
n_win(:,1) = floor(train_prctage * tot_n_win);  % Num. of train windows
n_win(:,2) = tot_n_win - n_win(:,1);  % Num of test windows

I = any(n_win<0, 2);
n_win(I,:) = [];
pnts(I) = [];
fnames(I) = [];

%% Adjust number of windows per epoch to match requested number
ndiff = sum(n_win) - [n_train n_test];
misswin = abs(diff(abs(ndiff)));
% we are missing windows: add one window from longest epochs
if misswin
    n_win(1:misswin, ndiff<0) = n_win(1:misswin, ndiff<0) + 1;
    ndiff = sum(n_win) - [n_train n_test];
end

if any(ndiff < 0)
    % Adjust number of train/test windows if not enough windows on one set
    k = abs(min(ndiff));
    [~, isort] = sort(sum(n_win, 2), 'descend');
    top_k_epochs = isort(1:k); % train
    n_win(top_k_epochs, ndiff<0) = n_win(top_k_epochs, ndiff<0) + 1;
    n_win(top_k_epochs, ndiff>0) = n_win(top_k_epochs, ndiff>0) - 1;
    ndiff = sum(n_win) - [n_train n_test];
end

if any(ndiff > 0)
    % If there is a surplus of k windows, substract one window from the top
    % k epochs with most number of windows
    [~, isort] = sort(n_win, 'descend');
    top_k_epochs = isort(1:ndiff(1),1); % train
    n_win(top_k_epochs,1) = n_win(top_k_epochs,1) - 1; 
    top_k_epochs = isort(1:ndiff(2),2); % test
    n_win(top_k_epochs,2) = n_win(top_k_epochs,2) - 1;     
end

I = all(n_win == [0 0], 2);
n_win(I,:) = [];
n_epochs = size(n_win, 1);
pnts(I) = [];
fnames(I) = [];
ndiff = sum(n_win) - [n_train n_test];
assert(all(ndiff==0), 'Something went wrong!')

%% Generate start indices for train and test windows on each epoch
if n_test == 0
    i_start = cell(n_epochs, 1);
else
    i_start = cell(n_epochs, 2);
end

for i_epoch = 1:n_epochs
    ni_train = n_win(i_epoch,1); % Num. of train windows in i-th epoch
    ni_test = n_win(i_epoch,2); ni = [ni_train ni_test];
    
    rand_start =  randperm(start_gap, 1);
    
    if overlap
        if n_test == 0
            i_start{i_epoch, 1} = uint32(rand_start + ...
                randperm(pnts(i_epoch)-start_gap, ni_train));
        else
            sublen = floor(pnts(i_epoch)/n_subdiv);
            offset = 0:sublen:sublen*(n_subdiv-1); % start indices of subintervals
            offset = offset(randperm(n_subdiv)); % shuffle
            offset = reshape(offset, n_subdiv/2, 2); % assign to train/test
            
            % Number of train/test windows on each subinterval
            sub_ni = floor(ni / (n_subdiv/2));
            sub_ni = repmat(sub_ni, n_subdiv/2, 1);
            kadd = ni - sum(sub_ni);
            sub_ni(1:kadd(1),1) = sub_ni(1:kadd(1),1) + 1;
            sub_ni(1:kadd(2),2) = sub_ni(1:kadd(2),2) + 1;
            
            i_start{i_epoch,1} = zeros(1, ni_train, 'uint32');
            i_start{i_epoch,2} = zeros(1, ni_test, 'uint32');
            win_cnt = [1 1];
            for i_sub = 1:n_subdiv/2
                
                % train
                ii_start = rand_start + randperm(sublen-start_gap, sub_ni(i_sub,1));
                ii_start = ii_start + offset(i_sub,1);
                i_start{i_epoch,1}(win_cnt(1):win_cnt(1)+sub_ni(i_sub,1)-1) = ...
                    uint32(ii_start);
                win_cnt(1) = win_cnt(1) + sub_ni(i_sub,1);
                
                % test
                ii_start = rand_start + randperm(sublen-start_gap, sub_ni(i_sub,2));
                ii_start = ii_start + offset(i_sub,2);
                i_start{i_epoch,2}(win_cnt(2):win_cnt(2)+sub_ni(i_sub,2)-1) = ...
                    uint32(ii_start);
                win_cnt(2) = win_cnt(2) + sub_ni(i_sub,2);
            end
        end
    else
        Lt = pnts(i_epoch) - start_gap; % Num. of points available for sampling
        max_step = floor(Lt/sum(ni));
        if max_step < winlen
            max_step = winlen;
        end        
        
        % Generate (ordered) start indexes that are at least winlen apart
        % (non-ovelapping windows)
        ii_start = rand_start:max_step:rand_start+max_step*(sum(ni)-1);                     
        ichoose = randperm(sum(ni)); % randomly shuffle ordered start indexes
        i_start{i_epoch,1} = uint32(ii_start(ichoose(1:ni(1))));  % train
        if n_test > 0
            i_start{i_epoch,2} = uint32(ii_start(ichoose(ni(1)+1:end))); % test
        end
    end
    
    % Make sure that we are not going out of boundaries
    max_i_start = max(horzcat(i_start{i_epoch,:}));
    assert(max_i_start+winlen-1 <= pnts(i_epoch), 'Something went wrong!');
    
    if n_test > 0
        % Make sure that there is no overlapping between train and test set
        rep_ni_test = repmat(i_start{i_epoch,2}', 1, ni_train);
        dist = abs(double(rep_ni_test)-double(i_start{i_epoch,1}));
        assert(all(dist(:) >= winlen), 'Train and test windows are overlapping!');
    end
        
end
end