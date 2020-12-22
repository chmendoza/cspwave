function dldata(data_dir, patient_id)
% Function to download data from the ieeg.org platform
%
% Parameters
% ----------
% data_dir (char)
%   Absolute path to folder where the EEG data set will be saved.
% patient (int)
%   Number of data set among a predefined list of data sets defined here.

if isdeployed
    patient_id = str2double(patient_id);
end

% data_dir = '/home/cmendoza/Research/sikmeans/LKini2019/data';
% patient_id = 22;
preictal_length = 3600 * 1e6; % usec away from EEC
horizon_lenght = 5 * 60 * 1e6; % usec away from EEC
interictal_gap = 4 * 3600 * 1e6; % usec away from ictal segment
interictal_length = 3600 * 1e6;
x_pat = '^x(?<id>\d+).mat$';

fprintf('Preferences directory: %s\n', prefdir);
max_samples_ieeg = 500 * 130 * 2000; % hz * channels * seconds
maxMemory = java.lang.Runtime.getRuntime.maxMemory; % bytes
freeAllocatedMemory = java.lang.Runtime.getRuntime.freeMemory;
usedMemory = java.lang.Runtime.getRuntime.totalMemory - ...
    freeAllocatedMemory;
freeMemory = maxMemory - usedMemory;
fprintf('Max Java heap space: %.3f MB\n', maxMemory/1e6);
fprintf('Free allocated Java heap space: %.3f MB\n', freeAllocatedMemory/1e6);
fprintf('Free reserved Java heap space: %.3f MB\n', freeMemory/1e6);


%% Labels to access datasets from ieee.org
dataset_names = {'HUP64_phaseII', 'HUP65_phaseII', 'HUP68_phaseII', ...
    'HUP70_phaseII', 'HUP73_phaseII', 'HUP74_phaseII', ...
    'HUP75_phaseII', 'HUP78_phaseII', 'HUP80_phaseII', ...
    'HUP82_phaseII', 'HUP83_phaseII', 'HUP86_phaseII', ...
    'HUP87_phaseII', 'HUP88_phaseII', 'HUP094_phaseII', ...
    'HUP105_phaseII', 'HUP106_phaseII', 'HUP107_phaseII', ...
    'HUP111_phaseII_D01', 'HUP111_phaseII_D02', ...
    'HUP116_phaseII', 'Study 012-2', 'Study 016', 'Study 017',...
    'Study 019', 'Study 022', 'Study 028', 'Study 029'};

%% Create session to access data from ieee.org
user = getenv('USER');
pwd_file = getenv('PWD_FILE');
session = IEEGSession(dataset_names{patient_id}, user, pwd_file);

%% Get from L. Kini et al., 2019: annotations, discarded electrodes, etc.
fname = getenv('LKini2019_metadata');
val = jsondecode(fileread(fname));
patients = rmfield(val.PATIENTS, 'Study020');
patient_labels = fieldnames(patients); 
patient_label = patient_labels{patient_id};

%% Remove ictal events with wrong annotations of EEC, UEO, or END
% EEC: Earliest electrographic change
% UEO: Unequivocal electrographic onset
% END: Termination of seizure
ictal = patients.(patient_label).Events.Ictal;
ictal = rmfield(ictal, 'x1000'); seiz_ids = fieldnames(ictal);
n_seiz = length(seiz_ids); seizTimes = zeros(n_seiz,3);
seizTimes(:,1) = structfun(@(x) x.SeizureEEC, ictal); % in seconds
seizTimes(:,2) = structfun(@(x) x.SeizureUEO, ictal);
seizTimes(:,3) = structfun(@(x) x.SeizureEnd, ictal);
I1 = any(seizTimes<0, 2);  I2 = any(diff(seizTimes,1,2)<0, 2); I = I1 | I2;
excluded_seiz = seiz_ids(I);
ictal = rmfield(ictal, excluded_seiz); 
seiz_ids(I) = []; seizTimes(I, :) = [];
[~, isort] = sort(seizTimes(:,1));
seizTimes = seizTimes(isort, :);
seizTimes = seizTimes * 1e6; % in usec: getvalues() units
seizTimes = seizTimes(:, [1 3]); % ictal segment: [EEC END]

ictal = structfun(@(x) rmfield(x, 'FILE'), ictal, 'UniformOutput', false);
I = structfun(@(x) isfield(x, 'STATUS'), ictal); nbad = sum(I);
seiz_str = repmat(' %s', 1, nbad); 
warn_msg = ['Seizures with STATUS field:', seiz_str, '\nDiscarding..\n'];
warning(warn_msg,  seiz_ids{I}); excluded_seiz = seiz_ids(I);
ictal = rmfield(ictal, excluded_seiz);
seiz_ids(I) = []; n_seiz = length(seiz_ids);

%% Reformat channel labels from ieeg.org
% to be consistent with labels from L. Kini et al., 2019.
% The latter are stored in patients, and have info about discarded electrodes.
match_expr = 'EEG (?<label>[a-zA-Z]*)(\s?)(0?)(?<number>\d*)-Ref';
replace_expr = '$<label>$<number>';
channel_labels = session.data.channelLabels(:,1);
channel_labels = regexprep(channel_labels, match_expr, replace_expr);

%% Get ids of channels, excluding ignored channels in LKini2019
n_channels = length(channel_labels);
channel_ids = 1:n_channels;
ignored_channels = patients.(patient_label).IGNORE_ELECTRODES;
ignored_ids = cellfun(@(x)find(strcmpi(channel_labels, x)), ignored_channels);
channel_ids = setxor(channel_ids, ignored_ids);
channel_labels = channel_labels(channel_ids);
n_channels = length(channel_labels);

%% Gather more metadata
srate = zeros(n_channels, 1); volt_factor = zeros(n_channels, 1);
n_samples = zeros(n_channels, 1);

for i_chan = 1:n_channels
    TimeSeriesDetails = session.data.rawChannels(channel_ids(i_chan)).get_tsdetails;
    srate(i_chan) = TimeSeriesDetails.getSampleRate;
    volt_factor(i_chan) = TimeSeriesDetails.getVoltageConversion;
    n_samples(i_chan) = TimeSeriesDetails.getNumberOfSamples;
end

%% Check that all the channels have the same sampling rate
assert(~any(abs(diff(srate))>0), "Channels don't have the same sampling rate");
srate = srate(1);
one_sample = 1e6/srate;

%% Check that all the channels have the same voltage conversion factor
assert(~any(abs(diff(volt_factor))>0), "Channels don't have the same voltage conversion factor");
% factor to multiply each sample by, to get the actual voltage reading (in uV).
volt_factor = volt_factor(1);

%% If not all the channels have the same number of samples, pick the min.
if any(abs(diff(n_samples))>0)
    n_samples = min(n_samples);
else
    n_samples = n_samples(1);
end

%% Build metadata struct
ignored_fields = {'RESECTION_IMAGE', 'Events', 'ELECTRODE_LABELS', 'IGNORE_ELECTRODES'};
metadata = rmfield(patients.(patient_label), ignored_fields);
metadata.ictal = ictal; metadata.channel_labels = channel_labels;
metadata.srate = srate; metadata.volt_factor = volt_factor;
metadata.n_samples = n_samples; metadata.snapName = session.data.snapName;

%% Make dirs
patient_dir = fullfile(data_dir, patient_label);
preictal_dir = fullfile(patient_dir, 'preictal');
interictal_dir = fullfile(patient_dir, 'interictal');
wd_dirs = {preictal_dir, interictal_dir};
if ~all(isfolder(wd_dirs))
    cellfun(@mkdir, wd_dirs)
else
    fprintf(['Skipping folder creation. Data already exists. ',...
        'Download will continue from its last check point...\n']);        
end

%% Create metadata file for patient
fname = fullfile(patient_dir, 'metadata.mat');
save(fname, 'metadata');

fprintf('Getting data from patient %s\n', patient_labels{patient_id})
fprintf('%d channels, %d seizures, %.2f samples/s\n', n_channels, n_seiz, srate);

%% Define maximum block size
% ieeg.org has a limit in the amount of data requested at a time, plus some
% java heap limits

max_samples_java_heap = floor(freeMemory * n_channels / 32); % double
max_block_size = floor(min(max_samples_ieeg, max_samples_java_heap)/n_channels);
max_block_size_us = floor(max_block_size/srate * 1e6); %in usec, per channel
fprintf('Max block size to be queried from ieeg.org: %.3f seconds\n', ...
    max_block_size_us/1e6);

%% Grab preictal segments

% Checkpoint: did we already downloaded some files?
file_list = dir(preictal_dir); file_list = {file_list.name};
x_files = regexp(file_list, x_pat, 'names');
x_files = x_files(~cellfun(@isempty, x_files));
x_id = cellfun(@(x)str2double(x.id), x_files);
if x_id
    cnt_file = max(x_id); %if stopped/failed, start from here
else
    cnt_file = 1;
end


t_start = [0; seizTimes(1:end-1,2)];
t_end = seizTimes(:,1);
inter_ictal_intervals = [t_start t_end];
% Some seizures are closer than preictal_zone_length+horizon_lenght
I = diff(inter_ictal_intervals, 1, 2) > (preictal_length+horizon_lenght);
preictal_intervals = [seizTimes(I,1)-preictal_length-horizon_lenght...
                      seizTimes(I,1)-horizon_lenght];
seiz_ids = seiz_ids(I); %attach seizure id to each preictal segment

n_preictal_segments = size(preictal_intervals, 1);
fprintf('Getting preictal segments...\n');
while cnt_file <= n_preictal_segments
    fname = sprintf('x%d.mat', cnt_file);
    fname = fullfile(preictal_dir, fname);
    mObj = matfile(fname, 'Writable', true);
    epoch_length = ceil(diff(preictal_intervals(cnt_file,:))*srate/1e6);    
    mObj.epoch = zeros(epoch_length, n_channels);
    mObj.seiz_id = seiz_ids{cnt_file};
    mObj.t_start = preictal_intervals(cnt_file,1); % usec
    mObj.t_end = preictal_intervals(cnt_file,2);
    
    n_blocks = ceil(epoch_length / max_block_size);
    t_start = preictal_intervals(cnt_file, 1); % in usec
    i_start = 1;    
    tStartRun = tic;
    for i_block = 1:n_blocks
        t_end = t_start + max_block_size_us;
        if t_end > preictal_intervals(cnt_file, 2)
            t_end = preictal_intervals(cnt_file, 2);
            block_size_us = t_end - t_start;            
        else
            block_size_us = max_block_size_us;
        end        
        aux = session.data.getvalues(t_start, block_size_us, channel_ids);
        block_size = size(aux, 1);
        mObj.epoch(i_start:i_start+block_size-1,:) = aux;
        i_start = i_start + block_size;
        t_start = t_end + one_sample;
    end
    mObj.epoch = mObj.epoch';
    tElapsed = toc(tStartRun);
    fprintf('Segment %d out of %d, in %.3f seconds\n', cnt_file, n_preictal_segments, tElapsed);
    cnt_file = cnt_file + 1;    
end

%% Grab interictal segments
tEnd = n_samples/srate*1e6; %final time in usec
t_start = [0; seizTimes(:,2)+interictal_gap];
t_end = [seizTimes(:,1)-interictal_gap; tEnd];
interictal_intervals = [t_start t_end];
% Some seizures are closer than 2*interictal_gap
I = diff(interictal_intervals,1,2) > 0;
interictal_intervals = interictal_intervals(I, :);

% Break into 1-hr segments because of RAM memory limitations
onehour_interictal_intervals = [];
n_interictal_segments = size(interictal_intervals, 1);
for ii = 1:n_interictal_segments
    t_start = interictal_intervals(ii,1):interictal_length+one_sample:interictal_intervals(ii,2);
    t_end = interictal_intervals(ii,1)+interictal_length:interictal_length+one_sample:interictal_intervals(ii,2);
    if t_end(end) < interictal_intervals(ii,2)
        t_end = [t_end interictal_intervals(ii,2)];
    end
    onehour_interictal_intervals = [onehour_interictal_intervals; [t_start' t_end']];
end
interictal_intervals = onehour_interictal_intervals;


% Checkpoint: did we already downloaded some files?
file_list = dir(interictal_dir); file_list = {file_list.name};
x_files = regexp(file_list, x_pat, 'names');
x_files = x_files(~cellfun(@isempty, x_files));
x_id = cellfun(@(x)str2double(x.id), x_files);
if x_id
    cnt_file = max(x_id); %if stopped/failed, start from here
else
    cnt_file = 1;
end

n_interictal_segments = size(interictal_intervals, 1);
fprintf('Getting interictal segments...\n');
while cnt_file <= n_interictal_segments
    fname = sprintf('x%d.mat', cnt_file);
    fname = fullfile(interictal_dir, fname);
    mObj = matfile(fname, 'Writable', true);
    epoch_length = ceil(diff(interictal_intervals(cnt_file,:))*srate/1e6);    
    mObj.epoch = zeros(epoch_length, n_channels);    
    mObj.t_start = interictal_intervals(cnt_file,1); % usec
    mObj.t_end = interictal_intervals(cnt_file,2);
    
    n_blocks = ceil(epoch_length / max_block_size);
    t_start = interictal_intervals(cnt_file, 1); % in usec
    i_start = 1;
    tStartRun = tic;
    for i_block = 1:n_blocks
        t_end = t_start + max_block_size_us;
        if t_end > interictal_intervals(cnt_file, 2)
            t_end = preictal_intervals(cnt_file, 2);
            block_size_us = t_end - t_start;
        else
            block_size_us = max_block_size_us;
        end
        aux = session.data.getvalues(t_start, block_size_us, channel_ids);
        block_size = size(aux, 1);
        mObj.epoch(i_start:i_start+block_size-1,:) = aux;
        i_start = i_start + block_size;
        t_start = t_end + one_sample;
    end
    mObj.epoch = mObj.epoch';
    tElapsed = toc(tStartRun);
    fprintf('Segment %d out of %d, in %.3f seconds\n', cnt_file, n_interictal_segments, tElapsed);
    cnt_file = cnt_file + 1;
end