% This function is to run the actual time frequency analysis based on the following paper:
% Meijer et al. (2016): Timing of beta oscillatory synchronization and temporal prediction of upcoming stimuli, NeuroImage
% 1. Load TFR data and motor performance data
% 2. Determine good subjects (With a minimum of 100 trials per condition)
% 3. cluster analysis for RP vs AP, SP vs AP and RS vs AS
% 4. correlation analysis (motor performance vs beta power for condition RP)
% 5. latency analysis: jackknife procedure following Miller, 1998

%% Set analysis variables
%necessary to run cluster analysis
load '\\192.168.32.1\fileshare\_projects\P3adaptiv\EEGcap_lay\lay_test.mat';

%subjects
sub={'03','07','08','09','11','12','13','14','15','19','20','21','22','23','25','26','27','28','29','30','34','35','36','37','38','39','43','44'};

%data
postfix='_meijer_2020_500';
pfad2='\\192.168.32.1\fileshare\_projects\Maren_Beta_TF\data\';
pfad='\\192.168.32.1\fileshare\_projects\Beta_P3adap\data\';

%conditions
con={'SI';'RI';'SR';'RR';'RA'};



%% -------------------------------------
% ------------ Load in data ------------
% --------------------------------------

clear('TFR_*');
i=1;

IBD_RI=[];

for o=1:length(sub)
    
    %load TFRs for all conditions
    for condi=1:length(con)
        
        filename=[pfad sprintf('S%s_%s%s.mat',sub{o},con{condi},postfix)];
        disp(filename)
        load(filename);
    end
    
    
    TFR_S=TFR_SI;
    TFR_S.powspctrm=TFR_SI.powspctrm-TFR_SR.powspctrm;
    TFR_S.subject=sub{o};
    TFRall_S{i}=TFR_S;
    
    TFR_R=TFR_RI;
    TFR_R.powspctrm=TFR_RI.powspctrm-TFR_RR.powspctrm;
    TFR_R.subject=sub{o};
    TFRall_R{i}=TFR_R;
    
    TFR_Ra=TFR_RA;
    TFR_Ra.powspctrm=TFR_RA.powspctrm-TFR_RR.powspctrm;
    TFR_Ra.subject=sub{o};
    TFRall_Ra{i}=TFR_Ra;
    
    TFR_SI.subject=sub{o};
    TFRall_SI{i}=TFR_SI;
    
    TFR_SR.subject=sub{o};
    TFRall_SR{i}=TFR_SR;
    
    TFR_RI.subject=sub{o};
    TFRall_RI{i}=TFR_RI;
    
    TFR_RR.subject=sub{o};
    TFRall_RR{i}=TFR_RR;
    
    TFR_RA.subject=sub{o};
    TFRall_RA{i}=TFR_RA;
    
    %load inter-beat interval deviations
    filename2=[pfad2 sprintf('S%s_tempomatching.mat',sub{o})];
    disp(filename2)
    load(filename2);
    IBD_RI=[IBD_RI tempoRARIRR(2)];
    
    i=i+1;
    
end

%% ---------------------------------------------------------------------------------------------
% ------------ Determine good subjects (with a minimum of 100 trials per condition) ------------
% ----------------------------------------------------------------------------------------------

lowsub=[];
goodsub=[];
threshold=100;


for j=1:length(sub)
    if (length(TFRall_RI{1,j}.trialinfo)<threshold) ||(length(TFRall_RR{1,j}.trialinfo)<threshold) ||(length(TFRall_SI{1,j}.trialinfo)<threshold) ||(length(TFRall_SR{1,j}.trialinfo)<threshold)
        disp(j);
        lowsub=[lowsub; TFRall_RR{1,j}.subject];
    else
        goodsub=[goodsub; j];
    end
end


%% 3.) Berechne Neighbors aus electroden positionen
load('\\192.168.32.1\fileshare\_projects\Beta_P3adap\Skripte2020\MarenElecPos.mat')
marenelec.chanpos=marenelec.elecpos;
ncfg=[];
ncfg.method='triangulation';
ncfg.elec          = marenelec;
neighbours       = ft_prepare_neighbours(ncfg);



goodsubj=goodsub;

%% 4.) calculate grand average
cfg=[];
cfg.keepindividual='yes';
cfg.parameter={'powspctrm'};
GA_SI=ft_freqgrandaverage(cfg,TFRall_SI{goodsubj});
GA_SR=ft_freqgrandaverage(cfg,TFRall_SR{goodsubj});
GA_RI=ft_freqgrandaverage(cfg,TFRall_RI{goodsubj});
GA_RR=ft_freqgrandaverage(cfg,TFRall_RR{goodsubj});
GA_RA=ft_freqgrandaverage(cfg,TFRall_RA{goodsubj});
GA_Ra=ft_freqgrandaverage(cfg,TFRall_Ra{goodsubj});
GA_S=ft_freqgrandaverage(cfg,TFRall_S{goodsubj});
GA_R=ft_freqgrandaverage(cfg,TFRall_R{goodsubj});
GA_RARR=GA_RA;
GA_RARR.powspctrm=GA_RA.powspctrm-GA_RR.powspctrm;

IBD=IBD_RI(goodsubj);



%% -----------------------------------------
% ------------ CLUSTER ANALYSIS ------------
% ------------------------------------------



cfg=[];
cfg.channel          = {'*'};
cfg.frequency        = [12 30];
cfg.avgoverfreq      = 'yes';
cfg.avgovertime      = 'no';
cfg.latency          = [-0.6   0];
cfg.avgoverchan      = 'no';
cfg.avgovertime      = 'no';
cfg.parameter        = 'powspctrm';
cfg.method           = 'montecarlo';
cfg.alpha            = 0.025;
cfg.numrandomization = 5000;
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.minnbchan        = 2;
cfg.neighbours       = neighbours;
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.ivar             = 2;
cfg.uvar             = 1;

cfg.design           = [1:size(GA_SI.powspctrm,1) 1:size(GA_SR.powspctrm,1); ones(1,size(GA_SI.powspctrm,1)) ones(1,size(GA_SR.powspctrm,1))*2 ];  %Designmatrix wird je nach Subjectanzahl automatisch angepasst

%% Rhythmic Pedaling versus Arrhythmic Pedaling (RP-AP)
statTallRIRR = ft_freqstatistics(cfg, GA_RI, GA_RR);
sig_timeRIRR=statTallRIRR.time(find(any(any(statTallRIRR.posclusterslabelmat==1,1),2)));
sig_timeRIRR_range=[sig_timeRIRR(1) sig_timeRIRR(end)];
sig_chanRIRR=statTallRIRR.label(find(any(any(statTallRIRR.posclusterslabelmat==1,3),2)));

%% Rhythmic Sitting versus Arrhythmic Sitting (RS-AS)
statTallSISR=ft_freqstatistics(cfg, GA_SI, GA_SR);
sig_timeSISR=statTallSISR.time(find(any(any(statTallSISR.posclusterslabelmat==1,1),2)));
sig_timeSISR_range=[sig_timeSISR(1) sig_timeSISR(end)];
sig_chanSISR=statTallSISR.label(find(any(any(statTallSISR.posclusterslabelmat==1,3),2)));

%% Rhythmic Pedaling versus Rhythmic Sitting (RP-RS)
cfg.latency=sig_timeRIRR_range;
statTallRISI=ft_freqstatistics(cfg, GA_RI, GA_SI);%
sig_timeRISI=statTallRISI.time(find(any(any(statTallRISI.posclusterslabelmat==1,1),2)));
sig_timeRISI_range=[sig_timeRISI(1) sig_timeRISI(end)];
sig_chanRISI=statTallRISI.label(find(any(any(statTallRISI.posclusterslabelmat==1,3),2)));

%% Arrhythmic Pedaling versus Arrhythmic Sitting (AP-AS) n.s.
statTallRRSR=ft_freqstatistics(cfg, GA_RR, GA_SR);%
sig_timeRRSR=statTallRRSR.time(find(any(any(statTallRRSR.posclusterslabelmat==1,1),2)));
sig_timeRRSR_range=[sig_timeRRRSR(1) sig_timeRRSR(end)];
sig_chanRRSR=statTallRRSR.label(find(any(any(statTallRRSR.posclusterslabelmat==1,3),2)));

%% Self-initiated Pedaling vs Arrhythmic Pedaling (SP vs AP)
cfg.latency=sig_timeRIRR_range;
statTallRARR=ft_freqstatistics(cfg, GA_RA, GA_RR);
sig_timeRARR=statTallRARR.time(find(any(any(statTallRARR.posclusterslabelmat==1,1),2)));
sig_timeRARR_range=[sig_timeRARR(1) sig_timeRARR(end)];
sig_chanRARR=statTallRARR.label(find(any(any(statTallRARR.posclusterslabelmat==1,3),2)));

%% Self-initiated Pedaling vs Rhythmic Pedaling (SP vs RP)
statTallRARI=ft_freqstatistics(cfg, GA_RA, GA_RI);
sig_timeRARI=statTallRARI.time(find(any(any(statTallRARI.posclusterslabelmat==1,1),2)));
sig_timeRARI_range=[sig_timeRARI(1) sig_timeRARI(end)];
sig_chanRARI=statTallRARI.label(find(any(any(statTallRARI.posclusterslabelmat==1,3),2)));


%% Self-initiated Pedaling versus  Rhythmic Sitting (SP-RS) n.s.
cfg.latency=sig_timeRARR_range;
cfg.channel=sig_chanRARR;
statTallRASI=ft_freqstatistics(cfg, GA_RA, GA_SI);



%% --------------------------------------------------------------------------------------------
% ------------ RP condition: Correlation motor performance/pre-stimulus beta power ------------
% ---------------------------------------------------------------------------------------------


cfg=[];
cfg.channel          = {'E17','E80','E100','E91','E119','E67','E71'}; %highest beta power RP-AP
cfg.avgoverchan      = 'no';
cfg.frequency        = [12 30];
cfg.avgoverfreq      ='no';
cfg.latency          = sig_timeRIRR_range;
cfg.avgovertime      = 'no';
cfg.parameter        = 'powspctrm';
cfg.method           = 'montecarlo';
cfg.alpha            = 0.05;
cfg.numrandomization = 5000;
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.correcttail      =  'alpha';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.minnbchan        = 2;
cfg.neighbours       = neighbours;
cfg.statistic        ='ft_statfun_correlationT';%
cfg.ivar             = 1;
cfg.design = [];
cfg.design(1,:) = IBD;
cfg.computestat    = 'yes';
freq_stat.correlation = ft_freqstatistics(cfg,GA_RI);


tmpstat = freq_stat.correlation;
tmpstat_timeRIRR = tmpstat.time(find(any(any(tmpstat.negclusterslabelmat==1,1),2)));
tmpstat_chanRIRR = tmpstat.label(find(any(any(tmpstat.negclusterslabelmat==1,3),2)));




