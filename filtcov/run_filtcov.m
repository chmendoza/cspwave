function run_filtcov(confpath)

% Script to run filtcov and save the covariance matrices to disk

if ~isdeployed
    addpath('../common');
    eeglab_path = getenv('EEGLAB_PATH');
    eeglab_functions = fullfile(eeglab_path, 'functions');
    firfilt_path = fullfile(eeglab_path, 'plugins/firfilt');    
    addpath(eeglab_path);
    addpath(genpath(eeglab_functions));
    addpath(genpath(firfilt_path));    
end

RunConf = jsondecode(fileread(confpath));

winlen = RunConf.winlen;
dirpath = RunConf.dirpath;
dfname = RunConf.dfname;
dpath = fullfile(dirpath, dfname);
mObj = matfile(dpath);
train_names = mObj.train_names;
train_indices = mObj.train_indices;

opts = {'dfnames', train_names,...
        'i_start', train_indices,...
        'winlen', winlen,...
        'cfname', RunConf.Sfname, ...
        'locutoff', RunConf.locutoff, ...
        'hicutoff', RunConf.hicutoff, ...
        'Fs', 512};

%% Train
filtcov(dirpath, opts{:});

end