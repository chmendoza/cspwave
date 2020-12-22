function prepdata(data_dir, patient)
% It calls prepdata_per_class() to preprocess the EEG data from a patient
%
% Parameters
% ----------
% data_dir (char)
%   Absolute path to folder with EEG data sets.
% patient (char)
%   Name of folder with preictal and interictal data from an EEG data set.
%   It has two subfolders, 'preictal' and 'interictal', with segments from 
%   each condition.

%% Delete files permanently
recycle('off');

%% Set up parpool (used by cleanLineNoise)
myCluster = parcluster('local');
myCluster.NumWorkers = str2double(getenv('SLURM_NTASKS_PER_NODE'));
myCluster.JobStorageLocation = getenv('TMPDIR');
mypool = parpool(myCluster, myCluster.NumWorkers);

if ~isdeployed
    addpath('../common');
    eeglab_path = getenv('EEGLAB_PATH');
    eeglab_functions = fullfile(eeglab_path, 'functions');
    firfilt_path = fullfile(eeglab_path, 'plugins/firfilt');
    prepline_path = fullfile(eeglab_path, 'plugins/PrepPipeline0.55.3');
    addpath(eeglab_path);
    addpath(genpath(eeglab_functions));
    addpath(genpath(firfilt_path));
    addpath(genpath(prepline_path));    
end

%% Load dataset metadata
patient_dir = fullfile(data_dir, patient);
fprintf('Getting data from patient %s\n', patient);
fname = fullfile(patient_dir, 'metadata.mat');
load(fname, 'metadata');

% multiply by this to get real voltage measure in uV
volt_factor = metadata.volt_factor; 

if isfield(metadata, 'orig_srate')
    metadata.srate = metadata.orig_srate;
end

n_channels = length(metadata.channel_labels);
srate = metadata.srate;
REF_Fs = 512; % Hz

%% Preprocessing paramaters
% cleanLineNoise parameters
lineNoiseIn = struct(...
    'lineNoiseChannels', 1:n_channels,...
    'Fs', srate, ...
    'lineFrequencies', [60, 120, 180],...
    'p', 0.01, ...
    'fScanBandWidth', 2, ...
    'taperBandWidth', 2, ...
    'taperWindowSize', 4, ...
    'taperWindowStep', 1, ...
    'tau', 100, ...
    'pad', 2, ...
    'fPassBand', [0 srate/2], ...
    'maximumIterations', 10);

min_epoch_length = 5; %seconds.
min_epoch_length = ceil(min_epoch_length * srate); % samples
buffer_length = 1; % minutes
buffer_length = ceil(buffer_length * 60 * srate); % samples

%% Set params struct

% One cycle of 40 Hz-rythm (low gamma). To detect dc chunks.
params.min_bad_chunk_length = floor(srate / 40);
% filtfilt imposes a limit in the minimum data length: 3*order+1.
% cleanLineNoise imposes a limit, depending on length of taper window and
% taper overlap.
params.min_epoch_length = min_epoch_length; % in samples
params.buffer_length = buffer_length; % in samples
% between 0 and 1, for rail-to-rail oscillation detection
params.freq_treshold = 0.05; % bigger is stricter
params.srate = srate;
params.locutoff = 1; % Hz;
params.lineNoiseIn = lineNoiseIn;
% To identify segments with 60-Hz hell, compare power around line freq.
% (line_band) with power in a contiguous band on the left (ref_band).
params.ref_band = [45 55]; % Hz
params.line_band = [55 65]; % Hz
params.REF_Fs = REF_Fs; % Hz;
% 70 uV. Ad hoc value, after visual inspection.
params.diff_threshold = 70/volt_factor;

%% Preictal dir
prepdata_per_class(patient_dir, 'preictal', params);

%% Preictal dir
prepdata_per_class(patient_dir, 'interictal', params);

if srate ~= REF_Fs
    %% Update metadata if resampling
    metadata.orig_srate = srate;
    metadata.srate = REF_Fs;
    save(fname, 'metadata');
end
delete(mypool);
end