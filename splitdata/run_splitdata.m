function run_splitdata(confpath)

RunConf = jsondecode(fileread(confpath));

data_dir = RunConf.data_dir;
winlen = RunConf.winlen; % window length
prctage = RunConf.prctage; % Percentage of data
train_prctage = RunConf.train_prctage; % Percentage of train data
overlap = RunConf.overlap; % True: Overlapping windows, False: non-overlapping
start_gap = RunConf.start_gap; % Gap for random start.
patient_label = RunConf.patient_label;
outfile = RunConf.outfile;  % Output file

conditions = {'preictal', 'interictal'};
n_cond = length(conditions);
for i_cond = 1:n_cond % each condition has its own folder
    
    fprintf('Splitting data in %s condition from patient %s\n',...
        conditions{i_cond}, patient_label);
        
    dirpath = fullfile(data_dir, patient_label, conditions{i_cond});
    
    [train_names, train_indices, test_names, test_indices] = ...
        splitdata(dirpath, winlen,...
                  'prctage', prctage,...
                  'train_prctage', train_prctage,...
                  'overlap', overlap,...
                  'start_gap', start_gap);
              
    outifile = fullfile(dirpath, outfile);
    save(outifile, 'train_names', 'train_indices', 'test_names', 'test_indices', '-v7.3');
    
end

end