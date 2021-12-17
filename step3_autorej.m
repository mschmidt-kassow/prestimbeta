% This function is to automatically remove artifacts:
% 1. Calculate channel dependencies (neighbours)
% 2. Load ICA cleaned data per subject
% 3. define bad trials and bad channels
% 4. repair bad channels
% 5. search for good trials in repaired channels
% 6. Save good trials

%% Set analysis variables

% channel neighbours
load('MarenElecPos.mat') %load electrode positions to determine channel neighbours
marenelec.chanpos=marenelec.elecpos;
ncfg=[];
ncfg.method='triangulation';
ncfg.elec          = marenelec;
ncfg.feedback = 'yes';
neighbours       = ft_prepare_neighbours(ncfg); %calculate channel neighbours

%subjects to read in
sub={'03','05','07','08','09','11','12','13','14','15','19','20','21','22','23','25','26','27','28','29','30','34','35','36','37','38','39','43','44'};


%variables
mvThresh=20; %Threshold for voltage differences between successive data points
mvThresh2=50; %Threshold for voltage differences between successive data points in repaired channels
maxArtifactsPerChannelPercent=60; %allowed percentage of artifacts per Channel that exceed the defined threshold

%% -------------------------------------
% ------------ Load in data ------------
% --------------------------------------

for o=1:length(sub)
    pfad='\\192.168.32.1\fileshare\_projects\Beta_P3adap\data\';
    loadname=[pfad sprintf('%s_clean_all.mat',sub{o})];
    eval(['load ' loadname ' ']);
    
    tmpdata=data;
    ntrl=length(tmpdata.trial);
    
  
    %Find maximum voltage differences for each channel
    level = zeros(length(tmpdata.label), ntrl);
    for j=1:ntrl
       level(:,j) = max(abs(diff(tmpdata.trial{j}')'),[],2);
    end
    BadTrialsPerChannel=sum(level>mvThresh,2)./(ones(size(level,1),1)*size(level,2)/100); %count differences that exceed predefined threhold
    BadChannel=find(BadTrialsPerChannel>maxArtifactsPerChannelPercent); %find channel(s) that has over 60% of bad trials
    
    %Repair "bad" channels via interpolation
    tmpcfg = [];
    tmpcfg.trials = 'all';
    tmpcfg.badchannel = tmpdata.label(BadChannel);
    tmpcfg.neighbours = neighbours;
    tmpcfg.elec = marenelec;
    repairdata = ft_channelrepair(tmpcfg, tmpdata);
    
    
    %Check repaired data
    ntrl=length(repairdata.trial);
    level2 = zeros(length(repairdata.label), ntrl);
    for j=1:ntrl
       level2(:,j) = max(abs(diff(repairdata.trial{j}')'),[],2); %count differences that exceed predefined threhold
    end
    
    GoodTrials=~any(level2>mvThresh2,1); %good trials are all trials that do not exceed the threshold
    
    
    %include only good trials
    cfg=[];
    cfg.trials=GoodTrials;
    data_clean = ft_preprocessing(cfg,repairdata);
           
    data=data_clean;
    
    % save clean data
    savename=[pfad sprintf('%s_autorej_0219.mat',sub{o})];
    eval(['save ' savename ' data -v7.3']);
end

%end of script. Got to step4_TFanalysis



