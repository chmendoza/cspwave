function prepdata_per_class(patient_dir, state, params)
% Preprocessing pipeline of EEG data
%
% prepdata_per_class(patient_dir, state, params) apply a preprocessing
% pipeline to data on patient_dir. The data is stored in a variable with
% name 'epoch', and has size [C T], with C and T being the number of EEG
% channels and time points, respectively. The variable 'epoch' is saved on
% files with pattern name 'x<epoch_id>.mat' (e.g., x13.mat). The epoch file
% also has the variables 'seiz_id' with the seizure id if the epoch is 
% preictal, and 't_start' with the start time of the epoch in the ieeg.org
% platform, in microseconds. After preprocessing each epoch file, new 
% preprocessed files are created with pattern name 'rx<chunk_id>.mat' (
% e.g., rx13.mat), and the raw epoch file is deleted. If state is 
% 'preictal', then the seizure id is saved into each chunk file. A new
% start time is also computed and saved for each chunk file.
%
% Parameters
% ----------
% patient_dir (char):  
%   Absolute path to folder with epoch data
% state (char):
%   If 'preictal', save seizure id that it is stored in on each chunk file
% params (struct):
%   Structure with fields
%       min_bad_chunk_length (int)
%           Minimum length of a chunk with constant value to be rejected
%       min_epoch_length (int)
%           filtfilt imposes a limit in the minimum data length that depends
%           on the order of the filter (3*order+1). cleanLineNoise also 
%           imposes a limit, depending on length of taper window and taper 
%           overlap.
%       buffer_length (int)
%           NUmber of samples to be added around segments with missing
%           values, same-amplitude chunks, and anomalous temporal spikes.
%       freq_treshold (float)
%           Threshold in relative frequency in [0 1] to discard segments
%           with rail-to-rail oscillation. bigger is stricter.
%       srate (float)
%           Sampling rate in Hz
%       locutoff (float)
%           Low-pass cut-off frequency in Hz
%       lineNoiseIn (struct)
%           Input structure for cleanLineNoise
%       ref_band (float)
%           Reference band to compare with spectral band power around line 
%           power. A vector of size [1 2] with values in Hz.
%       line_band (float)
%           Power line band. A vector of size [1 2] with values in Hz. A 
%           segment is rejected if its power in this band is higher to the 
%           power in ref_band. This deals with broad 60 Hz peak artifacts.
%       REF_Fs (float)
%           Reference sampling rate in Hz. The data is resampled to this 
%           value.
%       diff_threshold (float)
%           Threshold in absolute amplitude change between consecutive
%           samples to discard segments with anomalous temporal spkes.

x_pat = '^x(?<id>\d+).mat$';
rx_pat = '^rx(?<id>\d+).mat$';

[~, patient] = fileparts(patient_dir);
min_epoch_length = params.min_epoch_length;
buffer_length = params.buffer_length;
freq_treshold = params.freq_treshold;
srate = params.srate;
locutoff = params.locutoff;
lineNoiseIn = params.lineNoiseIn;
diff_threshold = params.diff_threshold;
ref_band = params.ref_band;
line_band = params.line_band;
REF_Fs = params.REF_Fs;
min_bad_chunk_length = params.min_bad_chunk_length;

one_sample = 1e6/srate; % in usec

wdir = fullfile(patient_dir, state);
file_list = dir(wdir); file_list = {file_list.name};

% Checkpoint: check if there are already processed chunks
rx_files = regexp(file_list, rx_pat, 'names');
rx_files = rx_files(~cellfun(@isempty, rx_files));
rx_id = cellfun(@(x)str2double(x.id), rx_files);
if rx_id
    cnt_file = max(rx_id) + 1; %if stopped/failed, start from here
else
    cnt_file = 1;
end

x_files = regexp(file_list, x_pat, 'names');
idx = ~cellfun(@isempty, x_files);
file_list = file_list(idx);
x_files = x_files(idx);
x_id = cellfun(@(x)str2double(x.id), x_files);
[~, isort] = sort(x_id);
file_list = file_list(isort);
n_files = length(file_list);

for i_file = 1:n_files
    
    tRunStart = tic;
    fprintf('Preprocessing %s\n', file_list{i_file});
    fpath = fullfile(wdir, file_list{i_file});
    mObj = matfile(fpath, 'Writable', true);
    [n_channels, n_samples] = size(mObj, 'epoch');
    
    %% Find NaN values (missing samples)
    nan_ind = find(any(isnan(mObj.epoch)));
    nan_ind = unique(nan_ind);
    nan_ind = find_chunk_indices(nan_ind, 1);
    nan_ind = interval_union(nan_ind);
    good_ind = get_good_segment_indices([1 n_samples], nan_ind);
    % Discard short segments
    chunk_length = diff(good_ind,1,2) + 1;
    good_ind = good_ind(chunk_length >= min_epoch_length, :);
    n_chunks = size(good_ind, 1);
    
    dc_ind = cell(n_chunks,1);
    for i_chunk = 1:n_chunks
        chunk_ind = good_ind(i_chunk,1):good_ind(i_chunk,2);
        %% Find dc segments (including all zeros). Apply to each channel indep.
        [channel, time_point] = find(diff(mObj.epoch(:,chunk_ind), 1, 2) == 0); %relative to chunk_ind
        time_point = chunk_ind(time_point)'; % Use absolute time reference
        channel = findgroups(channel);
        if ~isempty(channel)
            aux_dc_ind = splitapply(@(x){find_chunk_indices(x, min_bad_chunk_length)}, time_point, channel);
            dc_ind{i_chunk} = vertcat(aux_dc_ind{:});
        else
            dc_ind{i_chunk} = reshape([],[],2);
        end        
        dc_ind{i_chunk}(:,2) = dc_ind{i_chunk}(:,2) + 1; % diff leaves out the rightmost sample
        dc_ind{i_chunk} = add_buffer(dc_ind{i_chunk}, buffer_length, [1 n_samples]);        
    end
    dc_ind = cat(1, dc_ind{:});
    if isempty(dc_ind)
        dc_ind = reshape([],[],2);
    end
    
    %% Union of indices with NaN values and dc chunks
    bad_ind = union(nan_ind, dc_ind, 'rows');
    bad_ind = interval_union(bad_ind); % consolidate intersecting intervals
    
    %% Find indices of good segments
    good_ind = get_good_segment_indices([1 n_samples], bad_ind);
    
    %% Discard short segments
    chunk_length = diff(good_ind,1,2) + 1;
    good_ind = good_ind(chunk_length >= min_epoch_length, :);
    
    %% Discard segments with rail-to-rail oscillation
    n_chunks = size(good_ind, 1);
    fprintf(['Splitted into %d chunks after discarding segments with',...
             ' NaN values and long railed segments\n'], n_chunks);
    fprintf('Looking for chunks with rail-to-rail oscillation...\n');
    cnt_chunk = 1;
    israiled = @(x) x/sum(x) > freq_treshold;
    while cnt_chunk <= n_chunks
        chunk_ind = good_ind(cnt_chunk,1):good_ind(cnt_chunk,2);
        chunk_length = length(chunk_ind);
        t_start = mObj.t_start + (chunk_ind(1)-1) * one_sample;
        t_chunk_length = chunk_length * one_sample;
        
        G = kron(1:n_channels, ones(1, chunk_length));
        counts = splitapply(@(x){count_unique(x)}, reshape(mObj.epoch(:, chunk_ind)',1,[]),G);
        vals = splitapply(@(x){unique(x)}, reshape(mObj.epoch(:, chunk_ind)',1,[]),G);
        clear G;
        railed_chan_ind = cellfun(@(x)israiled(x), counts, 'UniformOutput', false);
        
        % Discard zero crossings
        nnz_ind = cellfun(@find, vals, 'UniformOutput', false);        
        railed_chan_ind = cellfun(@(x,y)x(y), railed_chan_ind, nnz_ind, 'UniformOutput', false);
        
        railed_chan_ind = cellfun(@any, railed_chan_ind);
        if any(railed_chan_ind) % bad chunk
            warn_msg = ['Railed chunk in patient %s, at %.3f seconds, ',...
                '%.3f seconds long. At channels', ...
                sprintf(' %d', find(railed_chan_ind))];
            warning(warn_msg, patient, t_start/1e6, t_chunk_length/1e6);
            good_ind(cnt_chunk,:) = [];
            n_chunks = n_chunks - 1;
            cnt_chunk = cnt_chunk - 1;
        end
        cnt_chunk = cnt_chunk + 1;
    end
    fprintf('Remaning chunks: %d\n', n_chunks);
    %% High-pass and line filtering, 60-Hz hell and anomalous peak removal
    cnt_chunk = 1;    
    while cnt_chunk <= n_chunks
        chunk_ind = good_ind(cnt_chunk,1):good_ind(cnt_chunk,2);
        chunk_length = length(chunk_ind);
        t_start = mObj.t_start + (chunk_ind(1)-1) * one_sample;
        t_chunk_length = chunk_length * one_sample;
        
        %% Get EEGLAB struct
        EEG = mat2EEGstruct(mObj.epoch(:, chunk_ind), srate);
        
        %% High pass filtering at 1 Hz (See Makoto's preprocessing pipeline)
        fprintf('High-pass filtering...\n');
        EEG = pop_eegfiltnew(EEG, 'locutoff', locutoff, 'plotfreqz', 0);
        
        %% Line filtering
        fprintf('Line filtering...\n');
        [EEG, ~] = cleanLineNoise(EEG, lineNoiseIn);
        
        %% Write to disk filtered data
        mObj.epoch(:, chunk_ind) = EEG.data;
        
        %% Discard chunks with 60-Hz hell
        fprintf('Looking for chunks with wide-band 60 Hz distortion...\n');
        ref_power = bandpower(EEG.data', srate, ref_band);
        line_power = bandpower(EEG.data', srate, line_band);
        has_wide_60Hz = line_power >= ref_power;
        if any(has_wide_60Hz)
            warn_msg = ['Chunk with wide-band 60 Hz distortion in ',...
                'patient %s, at %.3f seconds, %.3f seconds long. ',...
                'At channels', sprintf(' %d', find(has_wide_60Hz))];
            warning(warn_msg, patient, t_start/1e6, t_chunk_length/1e6);
            good_ind(cnt_chunk,:) = [];
            n_chunks = n_chunks - 1;
            continue;
        end
        fprintf('Remaning chunks: %d\n', n_chunks);
        %% Look for anomalous peaks        
        % i_sample are sample indexes within each chunk (EEG.data) for all channels.
        [chan, time_point] = find(abs(diff(EEG.data,1,2)) >= diff_threshold);
        time_point = unique(time_point);
        clear EEG;
        if time_point
            warn_msg = ['Chunk with anomalous peaks in ',...
                'patient %s, at %.3f seconds, %.3f seconds long. ',...
                'At channels', sprintf(' %d', unique(chan))];
            warning(warn_msg, patient, t_start/1e6, t_chunk_length/1e6);
            %% split further each chunk.
            % we need to convert those indexes within the whole epoch (mObj.epoch)
            time_point = chunk_ind(time_point)';
            % add buffer around sample with anomalous epoch
            time_point = repmat(time_point, 1, 2);
            bad_ind = add_buffer(time_point, buffer_length, [chunk_ind(1) chunk_ind(end)]);
            bad_ind = interval_union(bad_ind);
            good_subchunk_ind = get_good_segment_indices(...
                [chunk_ind(1) chunk_ind(end)], bad_ind);
            % Discard short subchunks
            subchunk_length = diff(good_subchunk_ind,1,2) + 1;
            good_subchunk_ind = good_subchunk_ind(subchunk_length >= min_epoch_length, :);
            if good_subchunk_ind
                %% Keep subchunks that are long enough
                n_subchunks = size(good_subchunk_ind, 1);
                good_ind = [good_ind(1:cnt_chunk-1,:);... % update indices
                    good_subchunk_ind;...
                    good_ind(cnt_chunk+1:end,:)];
                n_chunks = n_chunks - 1 + n_subchunks;
                cnt_chunk = cnt_chunk + n_subchunks - 1;
            end
        end
        %% Update counter
        cnt_chunk = cnt_chunk + 1;
    end
    fprintf('Final number of chunks: %d\n', n_chunks);
    %% Resample and save to disk each good chunk
    for i_chunks = 1:n_chunks
        chunk_ind = good_ind(i_chunks,1):good_ind(i_chunks,2);
        % Use 'rx' to avoid owerwritting 'x*.mat' files
        fname = sprintf('rx%d.mat', cnt_file);
        fname = fullfile(wdir, fname);
        chunk_mObj = matfile(fname, 'Writable', true);
        
        if srate ~= REF_Fs
            %% Resample from srate={500,499.9070} Hz to REF_Fs=512 Hz
            [P, Q] = rat(REF_Fs/srate);
            chunk_mObj.epoch = resample(mObj.epoch(:,chunk_ind)', P, Q)';
            one_sample = 1e6/REF_Fs;
        else
            chunk_mObj.epoch = mObj.epoch(:,chunk_ind);            
        end
        
        % Compute t_start and t_end for new chunk, in usec
        chunk_mObj.t_start = mObj.t_start + (chunk_ind(1)-1)*one_sample;
        chunk_mObj.t_end =  mObj.t_start + (chunk_ind(end)-1)*one_sample;
        if strcmpi(state, 'preictal')
            chunk_mObj.seiz_id  =  mObj.seiz_id;
        end
        cnt_file = cnt_file + 1;
    end    
    delete(fpath); % delete original file after chunks have been saved
    tElapsed = toc(tRunStart);
    fprintf('%s preprocessed and deleted after %.3f seconds\n', file_list{i_file}, tElapsed);
    fprintf('New file count: \n');
    disp(cnt_file);
end