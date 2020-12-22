function run_splitdata(confpath)

%% Set up parpool (used by cleanLineNoise)
myCluster = parcluster('local');
myCluster.NumWorkers = str2double(getenv('SLURM_NTASKS_PER_NODE'));
myCluster.JobStorageLocation = getenv('TMPDIR');
mypool = parpool(myCluster, myCluster.NumWorkers);

RunConf = jsondecode(fileread(confpath));

data_dir = '/lustre/scratch/cmendoza/sikmeans/LKini2019';

%mfname_frmt = 'monte%d_data_split_for_CSP_testing.mat';

conditions = {'preictal', 'interictal'};
n_cond = length(conditions);

winlen = RunConf.winlen; %512 * 2; % 2 seconds
n_train = RunConf.n_train; %100e3;
n_test = RunConf.n_test; % 0;
overlap = RunConf.overlap; % true;
n_subdiv = RunConf.n_subdiv; % 1;
goal = RunConf.goal;
n_monte = RunConf.n_monte;

patient_label = RunConf.patient_label;
parfor i_monte = 1:n_monte
    rng(i_monte);
    for i_cond = 1:n_cond % each condition has its own folder
        dirpath = fullfile(data_dir, patient_label, conditions{i_cond});
        [dfnames, i_start] = splitdata(dirpath, winlen, n_train(i_cond), ...
            n_test(i_cond),...
            overlap, n_subdiv);
        mfname = sprintf('monte%d_data_split_for_%s', i_monte, goal);
        fpath = fullfile(dirpath, mfname);
        mObj = matfile(fpath, 'Writable', true);
        mObj.fnames = dfnames;
        mObj.i_start = i_start;
    end
end
delete(mypool);
end