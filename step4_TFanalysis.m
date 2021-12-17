% This function is to run the actual time frequency analysis based on the following paper:
% Meijer et al. (2016): Timing of beta oscillatory synchronization and temporal prediction of upcoming stimuli, NeuroImage
% 1. Load cleaned data per subject
% 2. Resample data to 500 Hz
% 3. Check for remaining artifacts
% 4. Take only pre-stimulus interval
% 5. Padding of the data
% 6. TFR computation as described in Meijer et al, 2016
% 7. save TFR per condition

%% Set analysis variables

% subjects
subj={'03','05','07','08','09','11','12','13','14','15','19','20','21','22','23','25','26','27','28','29','30','34','35','36','37','38','39','43','44'};

mvThresh=80;

for s=1:length(subj)
    clear('TFR_*');
    
    %% -------------------------------------
    % ------------ Load in data ------------
    % --------------------------------------
    TFR_nrtrl=[];
    filename=['..\data\S' subj{s} '_autorej_0219.mat']; %my path to relevant input files
    disp(filename);
    load(filename);
    
    % reduce sampling rate to reduce processing costs
    
    cfg=[];
    cfg.detrend='yes';
    cfg.demean        = 'yes';
    cfg.baselinewindow  = 'all';
    cfg.trials = 'all';
    cfg.resamplefs      =  500 ;
    data = ft_resampledata(cfg, data)
    
    
    %Find maximum voltage differences for each channel
    ntrl=length(data.trial);
    level = zeros(length(data.label), ntrl);
    for j=1:ntrl
        level(:,j) = max(abs(diff(data.trial{j}')'),[],2);
    end
    
    
    % Reject trials with max > 80 µV
    [elecIDX,trialIDX] = find(level>mvThresh);
    trialIDX=unique(trialIDX);
    cfg = [];
    cfg.trials = 1:length(data.trial);
    cfg.trials(trialIDX)=[];
    data = ft_redefinetrial(cfg,data);
    
    % Take only pre-stimulus interval
    cfg                         = [];
    cfg.trials    = 'all';
    cfg.toilim    = [-1 0];
    data = ft_redefinetrial(cfg, data);
    
    %% preparation for TF analysis according to Meijer et al., 2016
    
    pfad='..\data\';
    fs = data.fsample;
    time_original =  data.time{1,1}; %Take copy of original timeline
    
    %Padding : Prior to TF analysis data were padded
    padtype = 'mirror'; %Apply mirror padding to data
    padlength = 2; %time (in seconds) up to which the trials are padded.
    ndatsample = size(data.trial{1,1},2);
    dattime = ndatsample/fs;
    npadsample = round((padlength-dattime)*fs);
    nprepadsample = floor(npadsample/2);
    npostpadsample = ceil(npadsample/2);
    prepadtimes = data.time{1,1}(1)+fliplr(-1*(1:nprepadsample))/fs;
    
    postpadtimes = data.time{1,1}(end)+(1:npostpadsample)/fs;
    
    for j=1:length(data.trial)
        data.trial{1,j} = ft_preproc_padding(data.trial{1,j}, padtype, nprepadsample, npostpadsample);
        data.time{1,j} = [prepadtimes data.time{1,j} postpadtimes];
    end
    
    %% ----------------------------------------
    % ------------ TFR Computation ------------
    % -----------------------------------------
    
    %1.) condition SI: corresponds to "RS" in the Ms (Rhythmic Sitting)
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.keeptrials   = 'yes';
    cfg.foi          = 6:(1/padlength):36;
    cfg.t_ftimwin    = 4./cfg.foi;
    cfg.toi          = time_original;
    
    cfg.trials=find(data.trialinfo(:,1)==12 & data.trialinfo(:,2)==4); %12= standard trials, 4= "condition 4"= Rhythmic Sitting
    TFR_SI = ft_freqanalysis(cfg, data);
    
    
    % Baseline normalization
    nr_of_samples = length(TFR_SI.time);
    trial_means = nanmean(TFR_SI.powspctrm,4);
    trial_means = repmat(trial_means,[1 1 1 nr_of_samples]);
    TFR_SI.powspctrm = ((TFR_SI.powspctrm - trial_means) ./ trial_means) * 100;          %percentage change relative to mean of trial (per frequency band)
    
    % average beta power over the trials
    TFR_SI.powspctrm = squeeze(mean(TFR_SI.powspctrm,1));
    TFR_SI.dimord = 'chan_freq_time';
    TFR_SI = rmfield(TFR_SI,'cumtapcnt');
    
    % save TFR for rhythmic sitting
    savename=[pfad sprintf('S%s_SI_meijer_2021_500.mat',subj{s})];
    eval(['save '  savename  ' TFR_SI* ' ]);
    
    clear('TFR_*');
    
    %2.) condition Sr: corresponds to "AS" in the Ms (arhythmic Sitting)
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.keeptrials   = 'yes';
    cfg.foi          = 6:(1/padlength):36;
    cfg.t_ftimwin    = 4./cfg.foi;
    cfg.toi          = time_original;
    
    cfg.trials=find(data.trialinfo(:,1)==12 & data.trialinfo(:,2)==5); %12= standard trials, 5= "condition 5"= arhythmic Sitting
    TFR_SR = ft_freqanalysis(cfg, data);
    
    nr_of_samples = length(TFR_SR.time);
    trial_means = nanmean(TFR_SR.powspctrm,4);
    trial_means = repmat(trial_means,[1 1 1 nr_of_samples]);
    TFR_SR.powspctrm = ((TFR_SR.powspctrm - trial_means) ./ trial_means) * 100;
    
    %Average over the trials
    TFR_SR.powspctrm = squeeze(mean(TFR_SR.powspctrm,1));
    TFR_SR.dimord = 'chan_freq_time';
    TFR_SR = rmfield(TFR_SR,'cumtapcnt');
    
    
    %save TFR for arhythmic Sitting
    savename=[pfad sprintf('S%s_SR_meijer_2021_500.mat',subj{s})];
    eval(['save '  savename  ' TFR_SR* ' ]);
    
    clear('TFR_*');
    
    %3.) condition RI: corresponds to "RP" in the Ms (rhythmic Pedaling)
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.keeptrials   = 'yes';
    cfg.foi          = 6:(1/padlength):36;
    cfg.t_ftimwin    = 4./cfg.foi;
    cfg.trials=find(data.trialinfo(:,1)==12 & data.trialinfo(:,2)==2); %12= standard trials, 2= "condition 2"= rhythmic Pedaling
    TFR_RI = ft_freqanalysis(cfg, data);
    
    
    nr_of_samples = length(TFR_RI.time);
    trial_means = nanmean(TFR_RI.powspctrm,4);
    trial_means = repmat(trial_means,[1 1 1 nr_of_samples]);
    TFR_RI.powspctrm = ((TFR_RI.powspctrm - trial_means) ./ trial_means) * 100;
    
    %Average over the trials
    TFR_RI.powspctrm = squeeze(mean(TFR_RI.powspctrm,1));
    TFR_RI.dimord = 'chan_freq_time';
    TFR_RI = rmfield(TFR_RI,'cumtapcnt');
    
    % save TFR for rhythmic Pedaling
    savename=[pfad sprintf('S%s_RI_meijer_2021_500.mat',subj{s})];
    eval(['save '  savename  ' TFR_RI* ' ]);
    
    clear('TFR_*');
    %4.) condition RR: corresponds to "AP" in the Ms (arhythmic Pedaling)
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.keeptrials   = 'yes';
    cfg.foi          = 6:(1/padlength):36;
    cfg.t_ftimwin    = 4./cfg.foi;
    cfg.toi          = time_original;
    cfg.trials=find(data.trialinfo(:,1)==12 & data.trialinfo(:,2)==3);%12= standard trials, 3= "condition 3"= arhythmic Pedaling
    TFR_RR = ft_freqanalysis(cfg, data);
    
    
    nr_of_samples = length(TFR_RR.time);
    trial_means = nanmean(TFR_RR.powspctrm,4);
    trial_means = repmat(trial_means,[1 1 1 nr_of_samples]);
    TFR_RR.powspctrm = ((TFR_RR.powspctrm - trial_means) ./ trial_means) * 100;          %percentage change relative to mean of trial (per frequency band)
    
    %Average over the trials
    TFR_RR.powspctrm = squeeze(mean(TFR_RR.powspctrm,1));
    TFR_RR.dimord = 'chan_freq_time';
    TFR_RR = rmfield(TFR_RR,'cumtapcnt');
    
    % save TFR for arhythmic Pedaling
    savename=[pfad sprintf('S%s_RR_meijer_2021_500.mat',subj{s})];
    eval(['save '  savename  ' TFR_RR* ' ]);
    clear('TFR_*');
    
    %5.) condition RA: corresponds to "SP" in the Ms (self-generated Pedaling)
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.keeptrials   = 'yes';
    cfg.foi          = 6:(1/padlength):36;
    cfg.t_ftimwin    = 4./cfg.foi;
    cfg.toi          = time_original;
    cfg.trials=find(data.trialinfo(:,1)==12 & data.trialinfo(:,2)==1); %12= standard trials, 1= "condition 1"= self-generated Pedaling
    TFR_RA = ft_freqanalysis(cfg, data);
    
    
    nr_of_samples = length(TFR_RA.time);
    trial_means = nanmean(TFR_RA.powspctrm,4);
    trial_means = repmat(trial_means,[1 1 1 nr_of_samples]);
    TFR_RA.powspctrm = ((TFR_RA.powspctrm - trial_means) ./ trial_means) * 100;          %percentage change relative to mean of trial (per frequency band)
    
    %Average over the trials
    TFR_RA.powspctrm = squeeze(mean(TFR_RA.powspctrm,1));
    TFR_RA.dimord = 'chan_freq_time';
    TFR_RA = rmfield(TFR_RA,'cumtapcnt');
    
    % save TFR for self-generated Pedaling
    savename=[pfad sprintf('S%s_RA_meijer_2021_500.mat',subj{s})];
    eval(['save '  savename  ' TFR_RA* ' ]);
    clear('TFR_*');
     
end

%end of script. Got to step5_clusteranalysis




