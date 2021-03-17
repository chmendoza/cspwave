function varargout = ...
    filtspatialfilt(W, dirpath, varargin)
% Spatially filter data extracted from data files
%
% Apply the spatial filters W to continous EEG segments. Each rx* file in
% dirpath is continous EEG segment that can be up to  1 hr long. Since each
% segment (file) might vary in length L, it returns the power (energy/L) of
% each CSP signal.
%
% If 'dfnames', 'i_start', and 'winlen' are provided, the function extracts
% windows of length winlen from files under dirpath and with names in
% dfnames. Apply the spatial filters W to those windows. It returns the
% energy (squared l2-norm) of the CSP signals.
%
% X = filtspatialfilt(W, dirpath, dfnames, i_start, winlen)
%   It returns all the spatiallly filtered windows, without bandpass
%   filtering. Rows are sample points.
%
% [X, dfnames, i_start] = filtspatialfilt(W, dirpath, dfnames, ...
%       i_start, winlen, 'choose', 'highenergy', 'k', k)
%   It returns the k spatially filtered windows with the highest energy per
%   spatial filter, the name of the epoch files, and their start indices on
%   the epoch variables.
%
% Parameters
% ----------
% W (double):
%   A matrix with size [C d], where C is the number of channels and d the
%   number of spatial filters.
% dirpath (str):
%   Absolute path of folder with data files. Each file has a variable
%   'epoch' of size [C T], with T being the number of time points.
%
% Optional name-value parameters
% ------------------------------
% 'dfnames' (str):
%   Cell array with the name of the files that will be used to extract the
%   windows. The pattern name is rx<epoch_id>.mat, with epoch_id being an
%   integer.
% 'i_start' (uint32):
%   Cell array with the start index of each window on each epoch file.
%   i_start{i}(j) is the start index of the j-th window in dfnames{i}.
%   size(i_start{i}) = [1 ni], where ni is the number of windows extracted
%   from dfnames{i}.
% 'winlen' (int):
%   Length of each window
% 'choose' (str)
%   Choose a subset of the windows. Choose all if not specified. Supported
%   options:
%     'highenergy': Choose the k filtered windows with the hightest energy
% 'k' (int)
%   The size of the subset of windows to be chosen
% 'locutoff' (float)
%   Low cut-off frequency. If not specified, don't do bandpass filtering.
%   of the bandpass filter, in Hz. Default: [].
% 'hicutoff' (float)
%   High cut-off frequency. If not specified, don't do bandpass filtering.
% 'Fs' (float)
%   Fs is the sampling rate in Hz. Needed for bandpass filtering.
%
% Optional output argument
% ------------------------
% X_csp (double):
%   A matrix of size [d T N], with d the number of spatial filters, T the number
%   of time points, and N the number of windows.
%   If not window-based processing, then it is a cell array of size [M 1],
%   with M the number of epochs.
%   It is only computed if the function is called with two output arguments:
%   [X, Px] = filtspatialfilt(W, dirpath).
% Ex (double):
%   A [c N] array with the squared l2-norm of the CSP signals
% Px (double):
%   The power of each CSP signal. A vector of size [d M].
% sub_dfnames (str)
%   Cell array with the name of the files that have the chosen subset of
%   windows. Size: [c 1]. size(sub_dfnames{i})=[1,M], with M being the
%   number of epoch files with the k chosen windows corresponding to the i-th
%   spatial filter. sub_dfnames{i}{j} is the name of the j-th epoch file
%   with windows chosen for the i-th spatial filter.
% sub_i_start (uint32)
%   Cell array with the start index of each chosen window on each epoch file.
%   sub_i_start{i}{j}(l) is the start index of the l-th window in
%   sub_dfnames{i}{j}. size(sub_i_start{i}{j}) = [1 ni], where ni is the
%   number of chosen windows in sub_dfnames{i}{j}.
% odfnames (char)
%   A cell array. odfnames{i}{j} is the file name from which window
%   X_csp(i,:,j) was taken.
% oi_start (uint32)
%   A cell array. oi_start{i}(j) is the start index of X_csp(i,:,j) inside
%   the epoch saved in file odfnames{i}{j}.
% winord (int)
%   A matrix of size [d k]. X_csp is NOT ordered according to the energy of
%   the CSP signals. windord(i,:) saves the order of the windows in X_csp
%   relative to their energy. For example, windord(i,13)=1 would mean that
%   X_csp(i,:,13) was the most energetic window after temporal and spatial
%   filtering.

%% Defaults for optional parameters
choose = []; k = []; locutoff = []; hicutoff = []; Fs = [];
dfnames = []; i_start = []; winlen = [];

%% Parse arguments
if nargin < 5
    error('Wrong number of arguments');
elseif nargin > 5
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
            case 'choose'
                choose = varargin{i_opt+1};
            case 'k'
                k = varargin{i_opt+1};
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


[n_chan, n_csp] = size(W);
keepXcsp = false;
if isempty(winlen)
    %% Process the entire epoch (rx* file)
    rx_pat = '^rx\d+.mat$';
    dfnames = dir(dirpath);
    dfnames = {dfnames.name};
    rx_files = regexp(dfnames, rx_pat);
    dfnames = dfnames(~cellfun(@isempty, rx_files));
    n_epochs = length(dfnames);
    
    Px = zeros(n_csp, n_epochs);
    if nargout == 2
        keepXcsp = true;
        X_csp = cell(n_epochs, 1);
    end
    i_start = cell(1, n_epochs);
else
    n_epochs = length(dfnames);
    if isempty(choose)
        X_csp = cell(n_epochs, 1); % [c T N]
    elseif strcmpi(choose, 'highenergy')
        if isempty(k)
            error('k must be specified');
        else
            Ex = cell(n_epochs, 1); % reduction in parfor
        end
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
max_pnts = max(pnts);
forStart = tic;

n_workers = str2double(getenv('SLURM_NTASKS_PER_NODE'));
if n_workers > n_epochs
    n_workers = n_epochs;
end
myCluster = parcluster('local');
myCluster.NumWorkers = n_workers;  %TODO: move this to the API
myCluster.JobStorageLocation = getenv('TMPDIR');
mypool = parpool(myCluster, myCluster.NumWorkers);


%% Initialize input variables
if isempty(winlen)
    if keepXcsp
        for i_epoch = 1:n_epochs
            X_csp{i_epoch} = ones(n_csp, pnts(i_epoch));
        end
    end
else
    n_winpep_ar = cellfun(@length, i_start);
    for i_epoch = 1:n_epochs
        if isempty(choose)
            X_csp{i_epoch} = ones(n_csp, winlen, n_winpep_ar(i_epoch));
        elseif strcmpi(choose, 'highenergy')
            Ex{i_epoch} = ones(n_csp, n_winpep_ar(i_epoch));
        end
    end
end

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
    
    %% Spatial filtering
    X = W'*X;
    
    if isempty(winlen)
        Px(:, i_epoch) = sum(X.^2,2)/n_samples;
        if keepXcsp
            X_csp{i_epoch} = X;
        end
    else
        %% Window-based processing
        for i_win = 1:n_winpep
            Xw = X(:, ii_start(i_win):ii_start(i_win)+winlen-1); % a window
            if isempty(choose)
                X_csp{i_epoch}(:,:,i_win) = Xw;
            elseif strcmpi(choose, 'highenergy')
                Ex{i_epoch}(:,i_win) = sum(Xw.^2,2);
            end
        end
    end
end
delete(mypool);
forEnd = toc(forStart);
fprintf('Main loop took %.3f seconds\n', forEnd);

if ~isempty(winlen)
    if isempty(choose)
        X_csp = cat(3, X_csp{:});
    elseif strcmpi(choose, 'highenergy')
        Ex = cat(2, Ex{:});
    end
end

if isempty(winlen)
    varargout{1} = Px;
    if keepXcsp
        varargout{2} = X_csp;
    end
else
    if isempty(choose)
        varargout{1} = X_csp;
    elseif strcmpi(choose, 'highenergy')
        X_csp = zeros(n_csp, winlen, k); %[c T k]
        sub_dfnames = cell(n_csp,1);
        sub_i_start = cell(n_csp,1);
        odfnames = cell(n_csp,1);        
        oi_start = cell(n_csp,1);
        winord = zeros(n_csp,k);
        ord = 1:k;
        G = kron(1:k, ones(1,n_epochs)); % to apply find(x,1) across columns of tf (see below)
        n_winpep = cellfun(@(x)length(x), i_start); % # of windows per file
        upbound = cumsum(n_winpep); % upper boundaries for i_win intervals. [n_epochs 1].
        prev_cum = [0; upbound(1:end-1)]'; % cnt_win at the start of each epoch iteration
        [~, isort] = sort(Ex, 2, 'descend');        
        for i_csp = 1:n_csp            
            % indices of the k windows with the highest energy in the i_csp-th
            % filter, respect to the whole set of windows:
            i_win = isort(i_csp, 1:k); % [1 k]
            tf = i_win <= upbound; % [n_epochs k]
            %index of dfnames{} for each chosen window
            i_file = splitapply(@(x)find(x,1),tf(:)',G);
            % index of each window whitin each epoch (file)
            ii_win = i_win - prev_cum(i_file);
            ui_file = unique(i_file);
            n_files = length(ui_file);
            sub_dfnames{i_csp} = dfnames(ui_file);
            sub_i_start{i_csp} = cell(1,n_files);
            cnt_win = 0;
            for ii_file = 1:n_files
                samefile = i_file==ui_file(ii_file);
                n_win = sum(samefile);                
                winord(i_csp,1+cnt_win:n_win+cnt_win)=ord(samefile);
                sub_i_start{i_csp}{ii_file} = ...
                    i_start{ui_file(ii_file)}(ii_win(samefile));
                cnt_win = cnt_win + n_win;
            end
            %% Call this function again, with some changes in varargin
            if i_csp == 1
                opts = varargin;
                opts = changeopts(opts, {'choose', 'k'}, {[],[]}); % Take all windows in recursive call
            end
            
            opts = changeopts(opts, 'dfnames', sub_dfnames{i_csp});
            opts = changeopts(opts, 'i_start', sub_i_start{i_csp}); 
            X_csp(i_csp, :, :) = filtspatialfilt(W(:,i_csp), dirpath, opts{:});
            odfnames{i_csp} = dfnames(i_file);
            oi_start{i_csp} = arrayfun(@(x,y) i_start{x}(y), i_file, ii_win);
        end
        varargout{1} = X_csp;
        varargout{2} = Ex;
        varargout{3} = sub_dfnames;
        varargout{4} = sub_i_start;
        varargout{5} = odfnames;
        varargout{6} = oi_start;
        varargout{7} = winord;
    end
end
end
