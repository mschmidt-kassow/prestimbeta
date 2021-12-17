% This function is to run the actual time frequency analysis based on the following paper:
% Meijer et al. (2016): Timing of beta oscillatory synchronization and temporal prediction of upcoming stimuli, NeuroImage
% 1. Load TFR data and motor performance data
% 2. Determine good subjects (With a minimum of 100 trials per condition)
% 3. Check for remaining artifacts
% 4. Take only pre-stimulus interval
% 5. Padding of the data
% 6. TFR computation as described in Meijer et al, 2016
% 7. save TFR per condition
% 8. cluster analysis for RP vs AP, SP vs AP and RS vs AS
% 9. correlation analysis (motor performance vs beta power for condition RP)
% 10. latency analysis: jackknife procedure following Miller, 1998

%% Set analysis variables
%necessary to run cluster analysis
load '\\192.168.32.1\fileshare\_projects\P3adaptiv\EEGcap_lay\lay_test.mat';

%subjects
sub={'03','05','07','08','09','11','12','13','14','15','19','20','21','22','23','25','26','27','28','29','30','34','35','36','37','38','39','43','44'};

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
load('\\192.168.32.1\fileshare\_projects\Beta_P3adap\Skripte0219\MarenElecPos.mat')
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


%% Self-initiated Pedaling vs Arrhythmic Pedaling (SP vs AP)
cfg.latency=sig_timeRIRR_range;
statTallRARR=ft_freqstatistics(cfg, GA_RA, GA_RR);
sig_timeRARR=statTallRARR.time(find(any(any(statTallRARR.posclusterslabelmat==1,1),2)));
sig_timeRARR_range=[sig_timeRARR(1) sig_timeRARR(end)];
sig_chanRARR=statTallRARR.label(find(any(any(statTallRARR.posclusterslabelmat==1,3),2)));









%% --------------------------------------------------------------------------------------------
% ------------ RP condition: Correlation motor performance/pre-stimulus beta power ------------
% ---------------------------------------------------------------------------------------------


cfg=[];
cfg.channel          = sig_chanRIRR;
cfg.avgoverchan      = 'no';
cfg.frequency        = [20 30];
cfg.avgoverfreq      ='yes';
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
timenegclus=find(any(any(tmpstat.negclusterslabelmat==1,1),2));


%% -------------------------------------------------------------------------------------
% ------------ Latency analysis: beta power peak; according to  Miller 1998 ------------
% --------------------------------------------------------------------------------------

%save relevant variables for jackknife procedure (latency analysis)
param_jackknife={};
param_jackknife.label={'RIRR';'SISR';'RARR'};
param_jackknife.time=[sig_timeRIRR_range;sig_timeSISR_range;sig_timeRARR_range];
param_jackknife.sigtimeRIRR=sig_timeRIRR;
param_jackknife.sigtimeSISR=sig_timeSISR;
param_jackknife.sigtimeRARR=sig_timeRARR;
param_jackknife.channelRIRR=[sig_chanRIRR];
param_jackknife.channelSISR=[sig_chanSISR];
param_jackknife.channelRARR=[sig_chanRARR];
save('param_jackknife_fixfreq.mat','param_jackknife');


%1.) Calculate grand average for each condition in significant time window
load('param_jackknife_fixfreq.mat');

cfg=[];
cfg.keepindividual='yes';
cfg.parameter={'powspctrm'};
cfg.toilim         = [param_jackknife.sigtimeSISR(1) param_jackknife.sigtimeSISR(end)];
GA_S=ft_freqgrandaverage(cfg,TFRall_S{goodsubj});
GA_SI=ft_freqgrandaverage(cfg,TFRall_SI{goodsubj});

cfg.toilim         = [param_jackknife.sigtimeRIRR(1) param_jackknife.sigtimeRIRR(end)];
GA_R=ft_freqgrandaverage(cfg,TFRall_R{goodsubj});
GA_RI=ft_freqgrandaverage(cfg,TFRall_RI{goodsubj});

cfg.toilim         = [param_jackknife.sigtimeRARR(1) param_jackknife.sigtimeRARR(end)] ;
GA_RR=ft_freqgrandaverage(cfg,TFRall_RR{goodsubj});
GA_RA=ft_freqgrandaverage(cfg,TFRall_RA{goodsubj});
GA_RaR=GA_RA;
GA_RaR.powspctrm=GA_RA.powspctrm-GA_RR.powspctrm;

%2.) Determine relevant electrodes, frequency and time per condition
%%% set PARAMETERs%%%%%
A=0.5; %Threshold

SISRall=eval(['GA_S.powspctrm']);
RIRRall=eval(['GA_R.powspctrm']);
RARRall=eval(['GA_RaR.powspctrm']);

nos=size(GA_S.powspctrm,1);

freqALL=find(GA_S.freq>=12 & GA_S.freq<30); %freqrange (same for all conditions)

%cluster channels for each condition
channelSISR=param_jackknife.channelSISR;
channelRIRR=param_jackknife.channelRIRR;
channelRARR=param_jackknife.channelRARR;

selecSISR=[]; %RS vs AS
for i=1:length(channelSISR)
    disp('****************');
    disp(channelSISR{i});
    disp('****************');
    a=find(strcmp(GA_SI.label,channelSISR{i}));
    selecSISR=[selecSISR, a];
end
selecRIRR=[]; %RP vs SP
for i=1:length(channelRIRR)
    disp('****************');
    disp(channelRIRR{i});
    disp('****************');
    a=find(strcmp(GA_SI.label,channelRIRR{i}));
    selecRIRR=[selecRIRR, a];
end
selecRARR=[]; %SP vs AP
for i=1:length(channelRARR)
    disp('****************');
    disp(channelRARR{i});
    disp('****************');
    a=find(strcmp(GA_SI.label,channelRARR{i}));
    selecRARR=[selecRARR, a];
end


%3.) Calculate time point of maximum power for cluster electrodes
%RS vs AS
[xSISR,y]=max(squeeze(nanmean(nanmean(nanmean(SISRall(:,selecSISR,freqALL,:),1),2),3)));
latSISR= GA_S.time(y)
%RP vs AP
[xRIRR,y]=max(squeeze(nanmean(nanmean(nanmean(RIRRall(:,selecRIRR,freqALL,:),1),2),3))); % die ist so spät
latRIRR= GA_R.time(y)

%SP vs AP
[xRARR,y]=max(squeeze(nanmean(nanmean(nanmean(RARRall(:,selecRARR,freqALL,:),1),2),3)));
latRARR= GA_RaR.time(y)

% 4.) Calculate Power Threshold
SchwelleRIRR= A*xRIRR; %RP vs AP
SchwelleSISR= A*xSISR;%RS vs AS
SchwelleRARR= A*xRARR; %SP vs AP

latri=[];%RP vs AP
for r = 1: length(sig_timeRIRR)
    if nanmean(nanmean(nanmean(RIRRall(:,bestRIRR,freqALL,r),2),3),1)> SchwelleRIRR == 1
        l=1;
    else
        l=0;
    end
    
    latri= [latri l];
end
min_Riall=GA_R.time(find(latri,1));

latsi=[]; %RS vs AS
for r = 1:length(sig_timeSISR)
    if nanmean(nanmean(nanmean(SISRall(:,bestSISR,freqALL,r),2),3),1)> SchwelleSISR == 1
        l=1;
    else
        l=0;
    end
    
    latsi= [latsi l];
end
min_Siall=GA_S.time(find(latsi,1));

latra=[]; %SP vs AP
for r = 1: length(sig_timeRARR)
    if nanmean(nanmean(nanmean(RARRall(:,bestRARR,freqALL,r),2),3),1)> SchwelleRARR == 1
        l=1;
    else
        l=0;
    end
    
    latra= [latra l];
end
min_Raall=GA_RaR.time(find(latra,1));


%6.) Actual Jackknife Proc: Average N-1
PowerSISR=[];
PowerRIRR=[];
PowerRARR=[];
Mins=[];

%Determin earliest timepoint where beta power is higher than threshold
for i=1:26;
    RIRRj=RIRRall;
    RARRj=RARRall;
    SISRj=SISRall;
     
    SISRj(i,:,:,:)=[];
    [x,y]=max(squeeze(nanmean(nanmean(nanmean(SISRj(:,selecSISR,freqALL,:),1),2),3)));
    y= GA_S.time(y);
    PowerSISR=[PowerSISR; x y];
    SchwelleSISR= A*x;
    
    RIRRj(i,:,:,:)=[];
    [x,y]=max(squeeze(nanmean(nanmean(nanmean(RIRRj(:,selecRIRR,freqALL,:),1),2),3)));
    y= GA_R.time(y);
    PowerRIRR=[PowerRIRR; x y];
    SchwelleRIRR= A*x; 
    
    RARRj(i,:,:,:)=[];
    [x,y]=max(squeeze(nanmean(nanmean(nanmean(RARRj(:,selecRARR,freqALL,:),1),2),3)));
    y= GA_RaR.time(y);
    PowerRARR=[PowerRARR; x y];
    SchwelleRARR= A*x; 
    
    ri=[];
    for r = 1:length(sig_timeRIRR)
        if nanmean(nanmean(nanmean(RIRRj(:,selecRIRR,freqALL,r),2),3),1)>= SchwelleRIRR == 1
            l=1;
        else
            l=0;
        end
        
        ri= [ri l];
    end
    min_Ri=GA_R.time(find(ri,1));
    
    ra=[];
    for r = 1:length(sig_timeRARR)
        if nanmean(nanmean(nanmean(RARRj(:,selecRARR,freqALL,r),2),3),1)>= SchwelleRARR == 1
            l=1;
        else
            l=0;
        end
        
        ra= [ra l];
    end
    min_Ra=GA_RaR.time(find(ra,1));
    
    si=[];
    for r = 1:length(sig_timeSISR)
        if nanmean(nanmean(nanmean(SISRj(:,selecSISR,freqALL,r),2),3),1)>= SchwelleSISR == 1
            l=1;
        else
            l=0;
        end
        
        si= [si l];
    end
    min_Si=GA_S.time(find(si,1));
    
    D_R=min_Ri-min_Si;
    D_Ra=min_Ra-min_Si;
    Mins = [Mins; min_Ri min_Si D_R min_Ra D_Ra];
end

%%7. statistics

nos=26;
% see Miller 1998: calculate standard error and then t-values from standard error
D_all=min_Riall-min_Siall; %Difference between 
DRa_all=min_Raall-min_Siall; %Difference between SP-AP and RS-AS

J=(sum(Mins(:,3)))/(nos); %RP-AP and RS-AS
JRa=(sum(Mins(:,5)))/(nos);%SP-AP and RS-AS

sumpow=sum(power((Mins(:,3)-J),2)); %RP-AP and RS-AS
sumpowRa=sum(power((Mins(:,5)-JRa),2));%SP-AP and RS-AS

Standarderror=sqrt(sumpow*(25/26)); %RP-AP and RS-AS
StandarderrorRa=sqrt(((sumpowRa*25)/26)); %SP-AP and RS-AS


Var=sumpow*((nos-1)/nos); %RP-AP and RS-AS
VarRa=sumpowRa*((nos-1)/nos); %SP-AP and RS-AS

Standardabweichung=sqrt(Var); %RP-AP and RS-AS
StandardabweichungRa=sqrt(VarRa); %SP-AP and RS-AS

% t-values
t_RiS=D_all/Standarderror %RP-AP and RS-AS
t_RaS=DRa_all/StandarderrorRa %SP-AP and RS-AS

alphalow = (0.05/3)/2;
alphaup = 1-(0.05/3)/2;

upp = tinv(alphaup,(nos-1));
low = tinv(alphalow,(nos-1));
CI_RiS=D_all+upp*Standarderror; % Confidence interval  RP-AP and RS-AS
CI_RaS=DRa_all+upp*StandarderrorRa; % confidence interval SP-AP and RS-AS
