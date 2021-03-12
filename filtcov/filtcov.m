function varargout = ...
    filtcov(dirpath, varargin)
% Compute covariance matrices of windows extracted from data files
%
% It extracts windows of length winlen from files under dirpath and with
% names in dfnames.
%
% filtcov(dirpath, fnames, i_start, winlen, 'cfname', cfname)
%   It computes the covariance matrices without bandpass filtering and save
%   them under dirpath with the name 'cfname.mat'.
%
% S = filtcov(dirpath, fnames, i_start, winlen)
%   It computes the covariance matrices without bandpass filtering and
%   return them as S, a CxCxN, with C the number of channels and N the
%   number of windows.
%
% Parameters
% ----------
% dirpath (str):
%   Absolute path of folder with data files. Each file has a variable
%   'epoch' of size [C T], with T being the number of time points.
%
% Optional name-value parameters
% -------------------
% dfnames (str):
%   Cell array with the name of the files that will be used to extract the
%   windows. The pattern name is rx<epoch_id>.mat, with epoch_id being an
%   integer.
% i_start (uint32):
%   Cell array with the start index of each window on each epoch file.
%   i_start{i}(j) is the start index of the j-th window in dfnames{i}.
%   size(i_start{i}) = [1 ni], where ni is the number of windows extracted
%   from dfnames{i}.
% winlen (int):
%   Length of each window
% 'cfname': [], or cfname (str)
%   If [], filtcov must be called with one output argument to return S. If
%   cfname, S will be saved under dirpath with name 'cfname.mat'. Default: [].
% 'locutoff': [], or locutoff (float)
%   If [], don't do bandpass filtering. locutoff is low cut-off frequency
%   of the bandpass filter, in Hz. Default: [].
% 'hicutoff': [], or hicutoff (float)
%   If [], don't do bandpass filtering. hicutoff is high cut-off frequency
%   of the bandpass filter, in Hz. Default: [].
% 'Fs': [], or Fs (float)
%   If [], don't do bandpass filtering. Fs is the sampling rate in Hz.
%   Default: [].
% 'ref': [], or 'CAR'
%   If [], don't do re-referencing. If 'CAR', do common average
%   re-referencing. Default: [].
% 'useGPU': true or false
%   If true, use GPU to do the main computations. Default: true.
%
% Optional output argument
% ------------------------
% S (double):
%   A matrix of size [C C N], with C the number of channels and N the
%   number of windows.

%% Defaults for optional parameters
dfnames = []; i_start = []; winlen = [];
cfname = []; locutoff = []; hicutoff = []; Fs = [];

%% Parse arguments
if nargin < 4
    error('Wrong number of arguments');
elseif nargin > 4
    n_varargin = length(varargin);
    assert(mod(n_varargin,2) == 0, 'Wrong number of optional arguments');
    for i_opt = 1:2:n_varargin
        switch varargin{i_opt}
            case 'dfnames'
                dfnames = varargin{i_opt+1};
            case 'i_start'
                i_start = varargin{i_opt+1};
            case 'winlen'
                winlen = varargin{i_opt+1};
            case 'cfname'
                cfname = varargin{i_opt+1};
            case 'locutoff'
                locutoff = varargin{i_opt+1};
            case 'hicutoff'
                hicutoff = varargin{i_opt+1};
            case 'Fs'
                Fs = varargin{i_opt+1};            
            otherwise
                error('Wrong optional arguments');
        end
    end
end

%% Get number of windows, channels and epochs
mObj = matfile(fullfile(dirpath, dfnames{1}));
n_chan = size(mObj, 'epoch', 1);
n_epochs = length(dfnames);

S = cell(1, n_epochs);
if isempty(winlen)
    %% Process the entire epoch (rx* file)
    rx_pat = '^rx\d+.mat$';
    dfnames = dir(dirpath);
    dfnames = {dfnames.name};
    rx_files = regexp(dfnames, rx_pat);
    dfnames = dfnames(~cellfun(@isempty, rx_files));
    n_epochs = length(dfnames);
    i_start = cell(1, n_epochs); % has to defined, o.w, parfor complains
else
    %% Compute variance over windows
    n_epochs = length(dfnames);    
    n_winpep_ar = cellfun(@length, i_start);
    for i_epoch = 1:n_epochs
        S{i_epoch} = zeros(n_chan,n_chan,n_winpep_ar(i_epoch));
    end
end

%% Configure parallel pool
pnts = zeros(n_epochs,1);
for i_epoch = 1:n_epochs
    %% Get length of each epoch
    fpath = fullfile(dirpath, dfnames{i_epoch});
    mObj = matfile(fpath);
    pnts(i_epoch) = size(mObj, 'epoch', 2);
end

forStart = tic;
n_workers = str2double(getenv('SLURM_NTASKS_PER_NODE'));
if n_workers > n_epochs
    n_workers = n_epochs;
end
myCluster = parcluster('local');
myCluster.NumWorkers = n_workers;
myCluster.JobStorageLocation = getenv('TMPDIR');
mypool = parpool(myCluster, myCluster.NumWorkers);


%% Main loop
parfor i_epoch = 1:n_epochs

    %% Load the data
    fpath = fullfile(dirpath, dfnames{i_epoch});
    mObj = matfile(fpath);
    X = mObj.epoch;
    fprintf('Epoch %d out of %d.\n', i_epoch, n_epochs);
    
    if isempty(winlen)
        [~, n_samples] = size(X);
    else
        ii_start = i_start{i_epoch};
        n_winpep = size(ii_start, 2); % windows per epoch
        fprintf('Processing %d windows in epoch %d.\n', n_winpep, i_epoch);
    end
    
    if ~isempty(Fs)
        %% EEGLAB structure init
        X = mat2EEGstruct(X, Fs);
        %% Band pass filtering
        X = pop_eegfiltnew(X, 'locutoff', locutoff, ...
            'hicutoff', hicutoff,'plotfreqz', 0);
        X = X.data;
    end
    X = X - mean(X, 2);
    
    if isempty(winlen)
        S{i_epoch} = X*X'/(n_samples-1);        
    else
        %% Window-based processing
        for i_win = 1:n_winpep
            Xw = X(:, ii_start(i_win):ii_start(i_win)+winlen-1); % a window
            S{i_epoch}(:,:,i_win) = Xw*Xw'/(winlen-1);
        end
    end
end
delete(mypool);
forEnd = toc(forStart);
fprintf('Main loop took %.3f seconds\n', forEnd);

S = cat(3, S{:});

if isempty(cfname)
    if nargout == 1
        varargout{1} = S;
    else
        error('Output argument not specified')
    end
else
    save(fullfile(dirpath, [cfname, '.mat']), 'S', '-v7.3');
end
end