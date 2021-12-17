%in this script we calculate motor performance data for each pedaling condition
%1.) Read in motor data and tone presentation data for each condition (SP, RP, AP)
%2.) Inter-Beat Deviation calculation following Leow et al, 2015

%subjects
sub={'03','05','07','08','09','11','12','13','14','15','19','20','21','22','23','25','26','27','28','29','30','34','35','36','37','38','39','43','44'};

condtext={'RA', 'RI', 'RR'}; %RA=SP in the Ms; RI=RP in the Ms; RR=AP in the Ms;

subIBD=[];

subRate=[];

for o=1:length(sub)
    i=1;
   
    pathname=['\\192.168.32.1\fileshare\_projects\P3adaptiv\data\Stud\S' sub{o} '\'];;
    filenames2={{ ['S' sub{o} '_Adap_Radeln_1.mat'], ['S' sub{o} '_Adap_Radeln_2.mat']},... %SP
        { ['S' sub{o} '_Iso_Radeln_1.mat'],  ['S' sub{o} '_Iso_Radeln_2.mat']},... %RP
        { ['S' sub{o} '_Rand_Radeln_1.mat'], ['S' sub{o} '_Rand_Radeln_2.mat']}}; %AP
    
    
   
    condiIBD=[];
    condiRate=[];   
    
    for condIDX=1:length(condtext)
        block=[];
        diffdiffall=[];
        rate=[];
        for blockIDX=1:2  % 2 blocks per condition
            fullfilename=[pathname filenames2{condIDX}{blockIDX}];
            disp(fullfilename);
            load(fullfilename);
            
            %motor
            TrittTime=data.TrittTime; %read in time points of motor behavior (pedaling)
            TrittDiff=diff(TrittTime); %calculate distance between successive rotations
            TrittRate=1./diff(TrittTime); %pedaling rate in Hz
            
            
            subject=filenames2{condIDX}{blockIDX};
            
            %tones
            TonTime=data.TonTime; %read in time points of acoustic presentation
            TonDiff=diff(TonTime);%calculate distance between successive tones
            TonRate=1./diff(TonTime); % acoustic presentation in Hz
          
            abTritt=min(find(TrittTime>TonTime(1))); %find time point of 1st rotation
            
            %determine analysis window: from first rotation + tone presentation until last tone + rotation
            
            if (length(TrittTime)<=length(TonTime))
                anaTritt=TrittTime(abTritt:end);
                anaTon=TonTime(1:length(anaTritt));
            elseif(length(TrittTime)>length(TonTime))
                
                anaTon=TonTime(1:end);
                anaTritt=TrittTime(abTritt:(length(anaTon)+(abTritt-1)));
            end
                        
            %Inter-Beat Deviation calculation following Leow et al, 2015
            IBD=mean(abs(diff(anaTritt(~isnan(anaTritt)))-diff(anaTon(~isnan(anaTon)))))/mean(diff(anaTon(~isnan(anaTon))));
            
            block=[block abs(IBD)]; %save IBD
            rate=[rate mean(TrittRate)]; %save speed
            i=i+1;
        end
        
        condi=mean(block);
        condiIBD=[condiIBD condi];
        Rate=mean(rate);
        condiRate=[condiRate Rate];
        
    end
    tempoRARIRR=[condiIBD];
    
    rateRARIRR=[condiRate];
    
    
    %save data to hard disk
    trittfile=sprintf('Z:\\_projects\\Maren_Beta_TF\\data\\%s_IBD.mat',sub{o});
    trittfile4=sprintf('Z:\\_projects\\Maren_Beta_TF\\data\\%s_speed.mat',sub{o});
        
    save(trittfile,'tempoRARIRR')
    save(trittfile4,'rateRARIRR')

end
%end of script.



