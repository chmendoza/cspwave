function run_csp(confpath)
% Function to run the main experiments in the paper. It calls filtspatialfilt
% to do the temporal and spatial filtering, and to retrieve the windows with
% highest energy.

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

n_csp = RunConf.n_csp;
winlen = RunConf.winlen;
method = RunConf.method;
dirpath = {RunConf.dircond1, RunConf.dircond2};
dfname = RunConf.dfname;
Sfname = RunConf.Sfname;
rfname = RunConf.rfname;
locutoff = RunConf.locutoff; % Hz
hicutoff = RunConf.hicutoff;

patient_dir = regexp(dirpath{1},filesep, 'split');
patient_label = patient_dir{end-1};
patient_dir = strjoin(patient_dir(1:end-1), filesep);

%% Absolute path to covariance matrices
S1fname = fullfile(dirpath{1}, Sfname);
S2fname = fullfile(dirpath{2}, Sfname);

X = cell(2, 2); % two conditions, train/test
Xtf = cell(2, 2);
Ex = cell(2, 2);
dfnames = cell(2, 2);
i_start = cell(2, 2);
odfnames = cell(2, 2);
oi_start = cell(2, 2);
winord = cell(2, 2);

opts = {'dfnames', [],...
    'i_start', [],...
    'winlen', [],...
    'choose', [],...
    'k', 15,...
    'locutoff', locutoff,...
    'hicutoff', hicutoff,...
    'Fs', []};

fprintf('Processing patient %s in band [%.3f %.3f] Hz\n', patient_label, ...
    locutoff, hicutoff);

%% Compute spatial filters
[W, y] = csp(S1fname, S2fname, 'k', n_csp, 'method', method);

% Output variables saved with splitdata()
splitvars = {{'train_names', 'train_indices'},...
             {'test_names', 'test_indices'}};

for i_cond =1:2
    for i_set = 1:2
        %% File names and start indices of train/test windows
        mObj = matfile(fullfile(dirpath{i_cond}, dfname));
        dfnames{i_cond,i_set} = mObj.(splitvars{i_set}{1}); % file names 
        i_start{i_cond,i_set} = mObj.(splitvars{i_set}{2}); % start indices
        
        %% Temporal filtering and window-based processing
        % For quantitative analysis: AUC and boxplot
        % choose k windows with highest energy
        opts = changeopts(opts, 'choose', 'highenergy');
        opts = changeopts(opts, 'Fs', 512);
        opts = changeopts(opts, 'winlen', winlen);
        opts = changeopts(opts, 'dfnames', dfnames{i_cond,i_set});
        opts = changeopts(opts, 'i_start', i_start{i_cond,i_set});
        [tXtf, tEx, tdfnames, ti_start, todfnames, toi_start, twinord] = filtspatialfilt(W, dirpath{i_cond}, opts{:});         
        Xtf{i_cond,i_set} = tXtf;
        Ex{i_cond,i_set} = tEx;
        dfnames{i_cond,i_set} = tdfnames;
        i_start{i_cond,i_set} = ti_start;
        odfnames{i_cond,i_set} = todfnames;
        oi_start{i_cond,i_set} = toi_start;
        winord{i_cond,i_set} = twinord;
        
        %% No temporal filtering and window-based
        % Run over the top k windows found in previous step
        % For qualitative analysis of waveforms
        opts = changeopts(opts, 'choose', []);
        opts = changeopts(opts, 'Fs', []);
        for i_csp = 1:n_csp
            opts = changeopts(opts, 'dfnames', dfnames{i_cond,i_set}{i_csp});
            opts = changeopts(opts, 'i_start', i_start{i_cond,i_set}{i_csp});
            X{i_cond,i_set}{i_csp} = filtspatialfilt(W(:,i_csp), dirpath{i_cond}, opts{:});
        end
        
        % Cat filtered windows along the csp dimension
        X{i_cond,i_set} = cat(1, X{i_cond,i_set}{:});
    end
end
resultspath = fullfile(patient_dir, rfname);
save(resultspath, 'W', 'y', 'Ex', 'Xtf', 'X', 'dfnames', 'odfnames', 'oi_start', 'winord', 'i_start', '-v7.3');
end