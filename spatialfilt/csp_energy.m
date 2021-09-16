function csp_energy(confpath)
% Compute all the CSP filters and the energy of all the spatially-filtered
% test data.

if ~isdeployed
    addpath('../common');
    addpath('../csp');
    eeglab_path = getenv('EEGLAB_PATH');
    eeglab_functions = fullfile(eeglab_path, 'functions');
    firfilt_path = fullfile(eeglab_path, 'plugins/firfilt');
    addpath(eeglab_path);
    addpath(genpath(eeglab_functions));
    addpath(genpath(firfilt_path));
end

RunConf = jsondecode(fileread(confpath));

i_csp = RunConf.i_csp;
seglen = RunConf.seglen;
winlen = RunConf.winlen;
method = RunConf.method;
dirpath = {RunConf.dircond1, RunConf.dircond2};
dfname = RunConf.dfname;
Sfname = RunConf.Sfname;
rfname = RunConf.rfname;
locutoff = RunConf.locutoff; % Hz
hicutoff = RunConf.hicutoff;

n_win_per_seg = seglen/winlen;
n_csp = length(i_csp);

patient_dir = regexp(dirpath{1},filesep, 'split');
patient_label = patient_dir{end-1};
patient_dir = strjoin(patient_dir(1:end-1), filesep);

%% Absolute path to covariance matrices
S1fname = fullfile(dirpath{1}, Sfname);
S2fname = fullfile(dirpath{2}, Sfname);

Ex = cell(2, 1); % two conditions, test

opts = {'dfnames', [],...
    'i_start', [],...
    'winlen', [],...    
    'locutoff', locutoff,...
    'hicutoff', hicutoff,...
    'Fs', []};

fprintf('Processing patient %s in band [%.3f %.3f] Hz\n', patient_label, ...
    locutoff, hicutoff);

%% Compute spatial filters
[W, ~] = csp(S1fname, S2fname, 'method', method);

for i_cond = 1:2
    
    %% File names and start indices of train/test windows
    mObj = matfile(fullfile(dirpath{i_cond}, dfname));
    dfnames = mObj.test_names;
    i_start = mObj.test_indices;
    
    %% Temporal filtering and window-based processing
    % For quantitative analysis: AUC and boxplot
    % choose k windows with highest energy    
    opts = changeopts(opts, 'Fs', 512);
    opts = changeopts(opts, 'winlen', seglen);
    opts = changeopts(opts, 'dfnames', dfnames);
    opts = changeopts(opts, 'i_start', i_start);
    X = filtspatialfilt(W(:,i_csp), dirpath{i_cond}, opts{:});
    
    n_seg = size(X, 3);
    n_win = n_win_per_seg * n_seg;
    X = permute(X, [2 3 1]); % (seglen, n_seg, n_csp)
    X = reshape(X, winlen, n_win, n_csp);
    Ex{i_cond} = squeeze(vecnorm(X, 2, 1)); % (n_win, n_csp)
    
end
resultspath = fullfile(patient_dir, rfname);
save(resultspath, 'W', 'Ex', '-v7.3');
end
