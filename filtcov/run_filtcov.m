function run_filtcov(confpath)

% Script to run filtcov and save the covariance matrices to disk
% It calls filtcov two times: one for train and one for test set

if ~isdeployed
    addpath('../common');
    eeglab_path = getenv('EEGLAB_PATH');
    eeglab_functions = fullfile(eeglab_path, 'functions');
    firfilt_path = fullfile(eeglab_path, 'plugins/firfilt');    
    addpath(eeglab_path);
    addpath(genpath(eeglab_functions));
    addpath(genpath(firfilt_path));    
end

% winlen = 512 * 2; % 2 seconds
RunConf = jsondecode(fileread(confpath));

winlen = RunConf.winlen;
dirpath = RunConf.dirpath;
dfname = RunConf.dfname;
dpath = fullfile(dirpath, dfname);
mObj = matfile(dpath);
fnames = mObj.fnames; 
i_start = mObj.i_start(:,1); % train

opts = {'dfnames', fnames,...
        'i_start', i_start,...
        'winlen', winlen,...
        'cfname', RunConf.Sfname, ...
        'locutoff', RunConf.locutoff, ...
        'hicutoff', RunConf.hicutoff, ...
        'Fs', 512};

%% Train
filtcov(dirpath, opts{:});

end