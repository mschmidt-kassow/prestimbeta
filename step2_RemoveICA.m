% This function is to remove ICA components related to blink activity:
% 1. Load preprocessed and appended EEG data per subject
% 2. Rereference EOG channel to get more pronounced blinks
% 3. Correlate ICA components with EOG channel
% 4. Reject components with r>0.3
% 5. Save clean data for subsequent analysis steps

%% Set analysis variables

thresh=0.3; %Threshold for correlation

%subjects to read in
sub={'03','05','07','08','09','11','12','13','14','15','19','20','21','22','23','25','26','27','28','29','30','34','35','36','37','38','39','43','44'};



%% -------------------------------------
% ------------ Load in data ------------
% --------------------------------------

for o=1:length(sub)
    cd '\\192.168.32.1\fileshare\_projects\Beta_P3adap\data\'; %my path to the preprocessed data
    loadname=sprintf('%s_comp_all.mat',sub{o});
    eval(['load ' loadname ' ']);
    
    
    
    %Rereference vEOG data
    cfg=[];
    cfg.channel={'E102'};
    cfg.reref={'E127'};
    dataeogreref = ft_preprocessing(cfg,dataeog);
    
    
    %Concatenate EOG trials to get the same file length (EOG and ICA components) for correlation
    cateog=[];
    smpPerTrial=size(comp.trial{1},2);
    catcomp=zeros(size(comp.trial{1},1),smpPerTrial*length(comp.trial));
    for i=1:length(comp.trial)
        tmp=dataeogreref.trial{i};
        cateog=[cateog tmp];
        tmp = comp.trial{i};
        catcomp(:,(i-1)*smpPerTrial+1:i*smpPerTrial)=tmp;
        disp(sprintf('concat data %d/%d',[i,length(comp.trial)]))
    end
    
    
    disp('starting corr');
    
    % Correlate component data with EOG data
    highcor=[];
    toreject=[];
    disp('corr...');
    cormatrix=corr(catcomp',cateog');
    
    thresh=mean(cormatrix)+3*std(cormatrix); %determine threshold
    
    highcor=find(abs(cormatrix)>=thresh); %check if threshold is above 0.3
    toreject=highcor; %if threshold >0.3 then mark as "to reject"
    
    %Reject marked components
    cfg = [];
    cfg.component = toreject;
    data = ft_rejectcomponent(cfg, comp);
    
    %Save clean data
    savename=sprintf('%s_clean_all.mat',sub{o});
    eval(['save ' savename ' data -v7.3']);
end

%end of script. Got to step3_autorej


