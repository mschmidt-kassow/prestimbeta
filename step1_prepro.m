% This function is to preprocess the EEG data:
% 1. Load raw EEG data
% 2. Get triggers
% 3. Make epochs with the given triggers
% 4. Automatic muscle artifact rejection
% 5. Run ICA
% 6. Save ICA components and EOG data (for subsequent correlation with ICA
% components)


%% Set analysis variables

%subjects to read in
sub={'03','05','07','08','09','11','12','13','14','15','19','20','21','22','23','25','26','27','28','29','30','34','35','36','37','38','39','43','44'};

%%%%%%%%%%%%%%experimental conditions%%%%%%%%%%%
%RA=Self-generated Pedaling (i.e. "SP" in the Ms)
%RI=Rhythmic Pedaling (i.e. "RP" in the Ms)
%RR=arrhythmic Pedalinf (i.e. "AP" in the Ms)
%SI=Rhythmic Sitting (i.e. "RS" in the Ms)
%SR=Arrhythmic Sitting (i.e. "AS" in the Ms)
condtext={'RA', 'RI', 'RR', 'SI', 'SR'};

%%%%to run the script, you need an electrode layout
load 'Z:\_projects\P3adaptiv\AuswertScripte\LayKonz35.mat';

prestim_ms=1000; %pre-stimulus interval of interest

%% To run script, you need the Field trip toolbox
% Add fieldtrip path
%restoredefaultpath;
%addpath('Z:\_projects\DPOAE\OAEASSR\scripts\fieldtrip-20120808');
%ft_defaults;


%% -------------------------------------
% ------------ Load in data ------------
% --------------------------------------
for o=1:length(sub)
    %Radeln RAND,Radeln Adap,Radel ISO,Sitzen ISO,Sitzen RAND
    pathname=['Z:\_projects\P3adaptiv\data\Stud\' sub{o} '\']; %my path
    %each condition consists of 2 files
    %Adap_Radeln=SP
    %Iso_Radeln=RP
    %Iso_Sitzen=RS
    %Rand_Sitzen=AS
    %Rand_Radeln=AP
    %!Please note that the order of the conditions has to match the order
    %in the "condtext"-Variable (line 22)
    
    filenames={{ [sub{o} '_Adap_Radeln_1.eeg'], [sub{o} '_Adap_Radeln_2.eeg']},...
        { [sub{o} '_Iso_Radeln_1.eeg'],  [sub{o} '_Iso_Radeln_2.eeg']},...
        { [sub{o} '_Rand_Radeln_1.eeg'], [sub{o} '_Rand_Radeln_2.eeg']},...
        {[sub{o} '_Iso_Sitzen_1.eeg'], [sub{o} '_Iso_Sitzen_2.eeg']},...
        {[sub{o} '_Rand_Sitzen_1.eeg'], [sub{o} '_Rand_Sitzen_2.eeg']},...
        { [sub{o} '_Iso_Radeln_1.eeg'],  [sub{o} '_Iso_Radeln_2.eeg']},...
        { [sub{o} '_Rand_Radeln_1.eeg'], [sub{o} '_Rand_Radeln_2.eeg']},...
        };
    
    
    anz_artifact_muscle=[]; %create variable for later muscle-artifacts
    alldata=[]; %create variable for appended data
    alldataeog=[]; % create variable for eog data
    
    
    datasetIDX=1; 
    %% --------------------------------------------------------------------
    % ------------ Append all condition files for each subject ------------
    % ---------------------------------------------------------------------
    for condIDX=1:length(condtext) %loop for each condition
        for blockIDX=1:length(filenames{condIDX}) %loop for each block per condition (N=2)
            fullfilename=[pathname filenames{condIDX}{blockIDX}];
            hdr = read_header(fullfilename); %load file
            
            cfg                         = [];
            cfg.dataset                 = fullfilename;
            cfg.trialdef.eventtype      = 'Stimulus';
            cfg.trialdef.prestim        = prestim_ms/1000; 
            cfg.trialdef.poststim       = 1; 
            trials_cfg = ft_definetrial(cfg); %read all Triggers (Events)
            
            %% --------------------------------------------
            % ------------ Get standard trials ------------
            % ---------------------------------------------
            
            trls=trials_cfg.trl;
            
            
            % take standard trials (12) only if there is no prior deviant (13)
            trl=trials_cfg.trl(trials_cfg.trl(:,4)>1,:);
            tmpidx=2:length(trl);
            trls=trl(tmpidx(find( trl(tmpidx,4)==12 & trl(tmpidx-1,4)~=13 )),:);
            
            
            trls=[trls  condIDX*ones(length(trls),1) ]; %Insert Column with condition label (1-5)
            %% --------------------------------------------
            % ------------ first preprocessing ------------
            % ---------------------------------------------
            cfg=[];
            cfg.trl = trls;
            cfg.channel    = {'E*','-E9'}; % electrode E9 has to be excluded due to technical issues
            cfg.dataset    = fullfilename;
            cfg.continuous = 'yes';
            
            cfg.reref      = 'yes'; %rereference to common average
            cfg.refchannel    = 'all';
            %
            
            cfg.demean     = 'yes';
            cfg.bpfilter = 'yes';
            cfg.bpfreq   = [0.1 150]; %band-pass filtering
            
            alldata{datasetIDX} = ft_preprocessing(cfg);
            
            %% -------------------------------------------------
            % ------------ muscle artifact detection------------
            % --------------------------------------------------
            cfg            = [];
            cfg.artfctdef.zvalue.channel = 'E*';
            cfg.artfctdef.zvalue.cutoff      = 20;  
            cfg.artfctdef.zvalue.trlpadding  = -0.2;
            cfg.artfctdef.zvalue.fltpadding  = 0.2;
            cfg.artfctdef.zvalue.artpadding  = 0;
            cfg.artfctdef.zvalue.bpfilter    = 'yes';
            cfg.artfctdef.zvalue.bpfreq      = [110 140];
            cfg.artfctdef.zvalue.bpfiltord   = 9;
            cfg.artfctdef.zvalue.bpfilttype  = 'but';
            cfg.artfctdef.zvalue.hilbert     = 'yes';
            cfg.artfctdef.zvalue.boxcar      = 0.2;
            cfg.artfctdef.zvalue.interactive = 'no';
            [cfg, artifact_muscle] = ft_artifact_zvalue(cfg,alldata{datasetIDX});
            
            %reject detected muscle artifacts
            cfg=[];
            cfg.artfctdef.reject = 'complete'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
            cfg.artfctdef.muscle.artifact = artifact_muscle; %
            alldata{datasetIDX} = ft_rejectartifact(cfg,alldata{datasetIDX});
            
            datasetIDX=datasetIDX+1; %count-up index
            anz_artifact_muscle=[anz_artifact_muscle size(artifact_muscle,1)]; %track number of muscle artifacts
            
        end
    end
    
    %append all datasets and blocks from each subject
    cfg=[];
    data=ft_appenddata(cfg,alldata{:});
    
    %extract vertical EOG channels to separate dataset
    cfg=[];
    cfg.channel    = {'E127','E102'};
    dataeog = ft_preprocessing(cfg,data);
    
    %run ICA without vEOG channels
    cfg = [];
    cfg.method = 'fastica';
    cfg.channel    = {'E*','-E127','-E102'}; %exclude vEOG channels
    comp = ft_componentanalysis(cfg,data);
    
    %save 1 file for each subject for further analysis
    savename=sprintf('../data/%s_comp_all.mat',sub{o});
    eval(['save ' savename ' dataeog comp -v7.3']);
    
    
end

%end of script. Got to step2_RemoveICA


