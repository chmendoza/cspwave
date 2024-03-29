function varargout = splitdata(dirpath, winlen, varargin)
% Pick multichannel windows of length winlen at random for training and testing
%
% [fnames, i_start] = splitdata(dirpath, winlen, n_train, n_test, varargin)
%
% Parameters
% ----------
% dirpath (str):
%   Absolute path to the data folder
% winlen (int):
%   Length of the instance (window)
%
% Optional parameters
% -------------------
% n_train (int):
%   Number of training instances
% n_test (int):
%   Number of test instances
% prctage (float):
%   Percentage of the data in the folder to be used. A number in [0 1].
% train_prctage (float):
%   Percentage of train data.  A number in [0 1].
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
% train_names (str):
%   cell array with the filenames of the epochs from where the train 
%   windows were picked.
% train_indices (uint32):
%   cell array with start index of the train windows within each epoch
%   (file). train_indices{i}(j) is the start index of the j-th window in 
%   the file train_names{i}.
% test_names (str):
%   cell array with the filenames of the epochs from where the test 
%   windows were picked.
% test_indices (uint32):
%   cell array with start index of the test windows within each epoch
%   (file). test_indices{i}(j) is the start index of the j-th window in 
%   the file test_names{i}.
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
n_train = []; n_test = []; prctage = 1; train_prctage = 0.8;
overlap = false; n_subdiv = 2; start_gap = winlen;

%% Parse arguments
if nargin < 2
    error('Wrong number of arguments');
elseif nargin > 2
    n_varargin = length(varargin);
    assert(mod(n_varargin,2) == 0, 'Wrong number of optional arguments');
    for i_opt = 1:2:n_varargin
        switch varargin{i_opt}
            case 'n_train'
                n_train = varargin{i_opt+1};
            case 'n_test'
                n_test = varargin{i_opt+1};
            case 'prctage'
                prctage = varargin{i_opt+1};
            case 'train_prctage'
                train_prctage = varargin{i_opt+1};
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

% Discard epochs that can not accomodate at least one window
ind = pnts < winlen+start_gap;

if any(ind)
    fprintf('Discarding %d out %d epochs that are not long enough\n', ...
        sum(ind), n_epochs);
    pnts(ind) = [];
    fnames(ind) = [];
    n_epochs = length(pnts);
end

if not(overlap)    
    %% ==== non-overlaping windows ============
    
    avail_pnts = pnts - start_gap;
    max_win = floor(avail_pnts/winlen);
    tot_max_win = sum(max_win);
    fprintf('Max. number of non-overlapping windows in %d epochs: %d\n',...
        n_epochs, tot_max_win);
    
    if ~isempty(n_train) && ~isempty(n_test)
        train_prctage = n_train / (n_train + n_test);
        test_prctage = 1 - train_prctage;
        tot_req_win = n_train + n_test;
        prctage = tot_req_win/tot_max_win;
        assert(tot_req_win <= tot_max_win, 'Not enough data!');
    else
        test_prctage = 1 - train_prctage;
        n_train = floor(train_prctage * prctage * tot_max_win);
        n_test = ceil(test_prctage * prctage * tot_max_win);
    end
    
    %% Generate max. num of window indices at random per epoch
    rand_start = randi(start_gap, n_epochs, 1);
    f = @(x,y) x + (1:winlen:y-winlen+1);
    win_ind = arrayfun(f, rand_start-1, avail_pnts, 'UniformOutput', false);
    win_pep = cellfun(@length, win_ind);
    
    %% Compute number of train windows to draw from each epoch
    win_train = floor(win_pep * prctage * train_prctage);
    rem_win = n_train - sum(win_train);
    % add 1 window to top epochs to achieve desired number of test windows
    win_train(1:rem_win) = win_train(1:rem_win) + 1;
    
    %% Compute number of test windows to draw from each epoch
    win_test = floor(win_pep * prctage * test_prctage);
    rem_win = n_test - sum(win_test);
    % epochs that have spare windows
    idx  = find(win_pep - win_train - win_test > 0);
    % add 1 window to top epochs to achieve desired number of test windows
    win_test(idx(1:rem_win)) = win_test(idx(1:rem_win)) + 1;    
    train_indices = cell(n_epochs, 1);
    if n_test > 0; test_indices = cell(n_epochs, 1); end
        
    for i_epoch = 1:n_epochs
        %% Pick windows at random
        
        % Pick indices for training windows
        idx = randperm(win_pep(i_epoch), win_train(i_epoch));
        train_indices{i_epoch} = win_ind{i_epoch}(idx);
        
        if n_test > 0
            % Discard train indices to ensure that test indices are different
            win_pep(i_epoch) = win_pep(i_epoch) - win_train(i_epoch);
            win_ind{i_epoch}(idx) = [];
            
            % Pick indices for test windows
            idx = randperm(win_pep(i_epoch), win_test(i_epoch));
            test_indices{i_epoch} = win_ind{i_epoch}(idx);
        end
        
    end
    
    % Remove epochs that did not contribute training windows
    idx = ~cellfun(@isempty, train_indices);    
    varargout{1} = fnames(idx);
    varargout{2} = cellfun(@uint32, train_indices(idx),'UniformOutput',false);
    fprintf(['Returning %d training non-overlapping windows chosen at ', ...
        'random from %d epochs\n'], n_train, sum(idx)); 
    
    if n_test > 0
        idx = ~cellfun(@isempty, test_indices);
        varargout{3} = fnames(idx);
        varargout{4} = cellfun(@uint32, test_indices(idx),'UniformOutput',false);
        fprintf(['Returning %d test non-overlapping windows chosen at ', ...
            'random from %d epochs\n'], n_test, sum(idx));
    end   
else
    %% ======== Overlapping windows  ================
    
    %% Computing max. number of windows
    % TODO: This is likely wrong, as it is not considering the case where
    % n_train > n_test: the subdivision for training should be longer.
    max_win = zeros(n_epochs, 1);
    %% epochs that can acommodate at least one window per sub-division.
    ind1 = pnts >= n_subdiv * (winlen+start_gap);
    divlen = pnts(ind1) - n_subdiv * start_gap;
    max_win(ind1) = divlen - winlen;
    
    %% epochs that can acommodate at least one window for two sub-divisions.
    ind2 = pnts >= 2 * (winlen+start_gap) & pnts < n_subdiv * (winlen+start_gap);
    divlen = pnts(ind2) - 2 * start_gap;
    max_win(ind2) = divlen - winlen;
    
    %% shorter epochs
    ind3 = pnts < 2 * (winlen+start_gap);
    divlen = pnts(ind3) - start_gap;
    max_win(ind3) = divlen - winlen;
    
    tot_max_win = sum(max_win);
    fprintf('Max. number of overlapping windows in %d epochs: %d\n',...
        n_epochs, tot_max_win);
    
    if ~isempty(n_train) && ~isempty(n_test)
        train_prctage = n_train / (n_train + n_test);
        test_prctage = 1 - train_prctage;
        tot_req_win = n_train + n_test;
        assert(tot_req_win <= tot_max_win, 'Not enough data!');
    else
        test_prctage = 1 - train_prctage;
        n_train = floor(train_prctage * prctage * tot_max_win);
        n_test = ceil(test_prctage * prctage * tot_max_win);
    end
    
    if n_test == 0
        avail_pnts = pnts - start_gap;
        rand_start = randi(start_gap, n_epochs, 1);
        f = @(x,y) x + (1:winlen:y-winlen+1);
        win_ind = arrayfun(f, rand_start-1, avail_pnts, ...
            'UniformOutput', false);
        win_pep = cellfun(@length, win_ind);
        cumwin = cumsum(win_pep);
        idx = cumwin <= n_train;
        idx = find(idx, 1, 'last');
        tot_win = cumwin(idx);
        if tot_win == n_train
            fnames = fnames(1:idx);
            i_start = win_ind(1:idx);
        else % n_train > tot_win
            rem_win = n_train - tot_win;
            fnames = fnames(1:idx+1);
            i_start = [win_ind(1:idx); {win_ind{idx+1}(1:rem_win)}];
        end
        return;
    end
    
    win_ind = cell(n_epochs, 2);
    for i_epoch = 1:n_epochs
        if any(i_epoch == find(ind1)) || any(i_epoch == find(ind2))
            if any(i_epoch == find(ind2))
                tn_subdiv = 1;
            else
                tn_subdiv = n_subdiv/2;
            end
            % epochs that can acommodate at least one window per sub-division.
            train_divlen = floor(train_prctage*pnts(i_epoch)/tn_subdiv);
            test_divlen = floor(test_prctage*pnts(i_epoch)/tn_subdiv);
            pair_divlen = train_divlen + test_divlen;
            pair_offset = 0:pair_divlen:pair_divlen*(tn_subdiv-1); % start indices of pair of subintervals
            sub_offsets = [0 train_divlen];
            shuffle_ind = randi(2, tn_subdiv); % Bernoulli with p=0.5, and values 1, and 2.
            train_offset = pair_offset + sub_offsets(shuffle_ind);
            test_offset = pair_offset + sub_offsets(3-shuffle_ind);
            win_div_ind = cell(tn_subdiv, 2);
            for i_div = 1:tn_subdiv-1
                rand_start = randperm(start_gap, 1);
                ii_start = rand_start + train_offset(i_div);
                win_div_ind{i_div, 1} = ii_start + randperm(train_divlen-start_gap-winlen+1);
                rand_start = randperm(start_gap, 1);
                ii_start = rand_start + test_offset(i_div);
                win_div_ind{i_div, 2} = ii_start + randperm(test_divlen-start_gap-winlen+1);
            end
            pair_divlen = pnts(i_epoch) - pair_divlen*(tn_subdiv-1);
            train_divlen = train_prctage * pair_divlen;
            test_divlen = pair_divlen - train_divlen;
            rand_start = randperm(start_gap, 1);
            ii_start = rand_start + train_offset(tn_subdiv);
            win_div_ind{tn_subdiv, 1} = ii_start + randperm(train_divlen-start_gap-winlen+1);
            rand_start = randperm(start_gap, 1);
            ii_start = rand_start + test_offset(tn_subdiv);
            win_div_ind{tn_subdiv, 2} = ii_start + randperm(test_divlen-start_gap-winlen+1);
            win_ind{i_epoch, 1} = horzcat(win_div_ind{:,1});
            win_ind{i_epoch, 2} = horzcat(win_div_ind{:,2});
        else
            % shorter epochs
            rand_start = randperm(start_gap, 1);
            win_ind{i_epoch, 1} = rand_start + randperm(pnts(i_epoch)-start_gap-winlen+1);
        end
    end
    nn = [n_train n_test];
    i_start = cell(n_epochs, 2);
    for i_set = 1:2 % train/test
        win_pep = cellfun(@length, win_ind(:,i_set));
        cumwin = cumsum(win_pep);
        idx = cumwin <= nn(i_set);
        idx = find(idx, 1, 'last');
        tot_win = cumwin(idx);
        if tot_win == nn(1)
            fnames = fnames(1:idx);
            i_start(1:idx, i_set) = win_ind(1:idx, i_set);
        else % n_train > tot_win
            rem_win = n_train - tot_win;
            fnames = fnames(1:idx+1);
            i_start(1:idx+1, i_set) = ...
                [win_ind(1:idx,i_set); {win_ind{idx+1,i_set}(1:rem_win)}];
        end
    end
end
end