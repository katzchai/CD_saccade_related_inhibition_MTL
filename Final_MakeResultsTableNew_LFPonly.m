function [ResultsTable]=F_MakeResultsTableNew_LFPOnly(Px,ResultsFolder,NameResulttable,Location,rawData,...
    electrodeTable,trialTable,CorrMacro,directionVector)
%% First see if neuron has adequate FR
tic
relevantSaccadeIndices = 1:height(trialTable);%saccade indices
relevantSaccadeIndices = relevantSaccadeIndices(~isinf(trialTable.SaccadeStartStop(:,1))...
    & ~isnan(trialTable.SaccadeStartStop(:,1)));%Saccade indices without blinks
NumberofSaccades=numel(relevantSaccadeIndices);%number of saccades
saccadeStartTimes = trialTable.SaccadeStartStop(relevantSaccadeIndices, 1);%saccade start times

%% Other specs of electrode/CorrMacro

% Find Specific electrode positions in the electrode Table

relevantElectrodeIndexMacro = find(electrodeTable.MacroMicro == 'MACRO'); %find Macro data on the electrode table

%% Directionality (get the directions of saccades and the number of directions)

directionVector = directionVector(~isinf(directionVector));

Left_sacc=find(directionVector == -1);nLeft=numel(Left_sacc);
Right_sacc=find(directionVector == 1);nRight=numel(Right_sacc);
%% For normalization: get the std of the LFP during experiment
raw_LFP=rawData.data(CorrMacro,:);%Get the LFP for that electrode for the entire experiment
filt_raw_LFP=eegfilt(raw_LFP,1000,0,50);%remove 60 Hz
%plot(filt_raw_LFP);% In case you want to plot
All_raw_lfp_macro_V=filt_raw_LFP./0.000000030517578125000001;%Convert to Hz according to Nlx (all experiments used the same range
% therefore this number did not vary.
std_raw_lfp_macro=std(All_raw_lfp_macro_V)%get std
%% LFP analysis
timeBeforeSaccade=2500;
timeafterSaccade=2500;


epochedDataLFPMacro = createEpochsWithDuration(rawData.data,...
    rawData.time, saccadeStartTimes, timeBeforeSaccade, timeafterSaccade, relevantElectrodeIndexMacro);

 
%Plot to evaluate LFP condition
x = -2500:1:2500;
lfpMacro=squeeze(epochedDataLFPMacro(:,CorrMacro,:));%epoched data

filt_lfp_macro=eegfilt(lfpMacro,1000,0,50);%I filtered the data 

%Change to Volts
raw_lfp_macro=filt_lfp_macro./0.000000030517578125000001;%converts to volts indicated in Nlx file

norm_lfp_macro=raw_lfp_macro./std_raw_lfp_macro;

mean_raw_LFPmacro=mean(raw_lfp_macro); %
mean_norm_LFPmacro=mean(norm_lfp_macro); %


%plot (mean(lfp),'k')
figure;
plot (x,mean_raw_LFPmacro,'k','LineWidth',2);
xlim([-500 500])
title (electrodeTable.ChannelName(relevantElectrodeIndexMacro(CorrMacro)));

%%%Provisional%% Test normailzation, usually keep commented
% 
% yyaxis left
% hold on
% plot (x,mean_raw_LFPmacro,'m','LineWidth',.5);
% 
% yyaxis right
% plot (x,mean_norm_LFPmacro,'g','LineWidth',2);
% box off

%% Evaluate if LFP is good or not
promptMacro = 'Is Macro good?  n=0; y=1;  ';
MacroERPgood = input(promptMacro);
if MacroERPgood == 0
StructFields_ifFail_LFP
    return
elseif MacroERPgood == 1
    promptFR = 'LFP is Ok '
end
%% Make all waveforms negative (some of the ERP's showed negative slopes and some positive)
% here I changed the polarity to keep all respones positive.
promptMacroPosNeg = 'waveform pos or neg? pos=0 neg=1;   ';
MacroERPPosNeg = input(promptMacroPosNeg)

if MacroERPPosNeg == 1
    mean_raw_LFPmacro=mean_raw_LFPmacro*-1; 
    mean_norm_LFPmacro=mean_norm_LFPmacro*-1; 
end

  
  %% RMS amplitude Note 
  resp_raw=mean_raw_LFPmacro(2510:3010);%We measured the RMS amplitude of the signal 10ms after saccade onset (to remove
  %artifact +500 ms
  resp_norm=mean_norm_LFPmacro(2510:3010);;
  
  rms_raw=rms(resp_raw);
  rms_norm=rms(resp_norm);%This value is the one we used for measurments and stats

  %% Directionality
  %% IPSI and CONTRA (define if it is 'contraversive' or 'ipsiversive' trials
ElectrodeName=electrodeTable.ChannelName(relevantElectrodeIndexMacro(CorrMacro));
Left="L"
isLeft = strfind(ElectrodeName,'L')

%Here I defined if the ERPs corresponded to ipsiversive or contraversive
%responses and separated it correspondingly. I am only going to comment the
%firs part. The second part is only a repetition of the other condition. 
%All measures are performed for raw and normalized signal although we only
%used the normalized for the paper. 
if isLeft == 1
ipsi_color = [0.3, 0.7, 0.9];
contra_color = [1 0 0];

raw_lfp_macro_Left=raw_lfp_macro(Left_sacc,:);
norm_lfp_macro_Left=norm_lfp_macro(Left_sacc,:);
raw_lfp_macro_Right=raw_lfp_macro(Right_sacc,:);
norm_lfp_macro_Right=norm_lfp_macro(Right_sacc,:);

mean_raw_LFPmacro_Left=mean(raw_lfp_macro_Left); 
mean_norm_LFPmacro_Left=mean(norm_lfp_macro_Left); 
mean_raw_LFPmacro_Right=mean(raw_lfp_macro_Right); 
mean_norm_LFPmacro_Right=mean(norm_lfp_macro_Right); 

sacc_start=0;

%plot (mean(lfp),'k')
figure;
plot (x,mean_raw_LFPmacro_Left,'Color',ipsi_color,'LineWidth',2);
hold on
plot (x,mean_raw_LFPmacro_Right,'Color',contra_color,'LineWidth',2);

xlim([-500 500])
title (electrodeTable.ChannelName(relevantElectrodeIndexMacro(CorrMacro)));

%% Make all waveforms negative (again here I converted all ERP's to positive polarity)
promptMacroPosNeg = 'waveform pos or neg? pos=0 neg=1;   ';
MacroERPPosNeg = input(promptMacroPosNeg)

if MacroERPPosNeg == 1
    mean_raw_LFPmacro_Left=mean_raw_LFPmacro_Left*-1; 
    mean_norm_LFPmacro_Left=mean_norm_LFPmacro_Left*-1; 
    mean_raw_LFPmacro_Right=mean_raw_LFPmacro_Right*-1; 
    mean_norm_LFPmacro_Right=mean_norm_LFPmacro_Right*-1; 
end

%%

Direction_ResultsFolder= 'D:\Matlab\neurol\data\Saccade_New\RsultsSaccade_v3\';
FolderNameDirection_LFP= '\Figures\Direction\LFP';

GenAdj_LFP= "ERP_Direction";
Pxfolder=[cell2mat(Px)];

 cd([Direction_ResultsFolder Pxfolder FolderNameDirection_LFP])
    format_fig=".fig";
    format_jpg=".jpg";
    format_eps=".eps";
    space="_";
    figname_fig=strcat(num2str(GenAdj_LFP),space,Px,space, electrodeTable.ChannelName(relevantElectrodeIndexMacro(CorrMacro)), format_fig);
    figname_jpg=strcat(num2str(GenAdj_LFP),space,Px,space, electrodeTable.ChannelName(relevantElectrodeIndexMacro(CorrMacro)), format_jpg);
    figname_eps=strcat(num2str(GenAdj_LFP),space,Px,space, electrodeTable.ChannelName(relevantElectrodeIndexMacro(CorrMacro)), format_eps);
    
    saveas(gcf,figname_fig);
    saveas(gcf,figname_jpg);
    saveas(gcf,figname_eps);
    close all
   %% Get general rms amplitude

  resp_raw_Left=mean_raw_LFPmacro_Left(2510:3010);
  resp_norm_Left=mean_norm_LFPmacro_Left(2510:3010);
  
  rms_raw_Left=rms(resp_raw_Left);
  rms_norm_Left=rms(resp_norm_Left);

  resp_raw_Right=mean_raw_LFPmacro_Right(2510:3010);
  resp_norm_Right=mean_norm_LFPmacro_Right(2510:3010);
  
  rms_raw_Right=rms(resp_raw_Right);
  rms_norm_Right=rms(resp_norm_Right);
  %%
  
  ipsilateral_side="Left";
  contralateral_side="Right";

 % Rename %To do analysis later merging all the ipsilateral responses and
 % all the contralateral responses.
    mean_raw_LFPmacro_Ipsi=mean_raw_LFPmacro_Left; 
    mean_norm_LFPmacro_Ipsi=mean_norm_LFPmacro_Left; 
    mean_raw_LFPmacro_Contra=mean_raw_LFPmacro_Right; 
    mean_norm_LFPmacro_Contra=mean_norm_LFPmacro_Right; 
     
 
    resp_raw_Ipsi=resp_raw_Left;
    resp_norm_Ipsi=resp_norm_Left;
    rms_raw_Ipsi=rms_raw_Left;
    rms_norm_Ipsi=rms_norm_Left;
    resp_raw_Contra=resp_raw_Right;
    resp_norm_Contra=resp_norm_Right;
    rms_raw_Contra=rms_raw_Right;
    rms_norm_Contra=rms_norm_Right;
elseif isempty (isLeft) %% Here you do the same in the case that it is an electrode from the right hemispere
ipsi_color = [0.3, 0.7, 0.9];
contra_color = [1 0 0];

raw_lfp_macro_Right=raw_lfp_macro(Right_sacc,:);
norm_lfp_macro_Right=norm_lfp_macro(Right_sacc,:);
raw_lfp_macro_Left=raw_lfp_macro(Left_sacc,:);
norm_lfp_macro_Left=norm_lfp_macro(Left_sacc,:);

mean_raw_LFPmacro_Right=mean(raw_lfp_macro_Right); 
mean_norm_LFPmacro_Right=mean(norm_lfp_macro_Right); 
mean_raw_LFPmacro_Left=mean(raw_lfp_macro_Left); 
mean_norm_LFPmacro_Left=mean(norm_lfp_macro_Left); 

sacc_start=0;

%plot (mean(lfp),'k')
figure;
plot (x,mean_raw_LFPmacro_Right,'Color',ipsi_color,'LineWidth',2);
hold on
plot (x,mean_raw_LFPmacro_Left,'Color',contra_color,'LineWidth',2);

xlim([-500 500])
title (electrodeTable.ChannelName(relevantElectrodeIndexMacro(CorrMacro)));

%% Make all waveforms negative
promptMacroPosNeg = 'waveform pos or neg? pos=0 neg=1;   ';
MacroERPPosNeg = input(promptMacroPosNeg)

if MacroERPPosNeg == 1
    mean_raw_LFPmacro_Right=mean_raw_LFPmacro_Right*-1; %4struct
    mean_norm_LFPmacro_Right=mean_norm_LFPmacro_Right*-1; %4struct
    mean_raw_LFPmacro_Left=mean_raw_LFPmacro_Left*-1; %4struct
    mean_norm_LFPmacro_Left=mean_norm_LFPmacro_Left*-1; %4struct
end

%%

Direction_ResultsFolder= 'D:\Matlab\neurol\data\Saccade_New\RsultsSaccade_v3\';
FolderNameDirection_LFP= '\Figures\Direction\LFP';

GenAdj_LFP= "ERP_Direction";
Pxfolder=[cell2mat(Px)];

 cd([Direction_ResultsFolder Pxfolder FolderNameDirection_LFP])
    format_fig=".fig";
    format_jpg=".jpg";
    format_eps=".eps";
    space="_";
    figname_fig=strcat(num2str(GenAdj_LFP),space,Px,space, electrodeTable.ChannelName(relevantElectrodeIndexMacro(CorrMacro)), format_fig);
    figname_jpg=strcat(num2str(GenAdj_LFP),space,Px,space, electrodeTable.ChannelName(relevantElectrodeIndexMacro(CorrMacro)), format_jpg);
    figname_eps=strcat(num2str(GenAdj_LFP),space,Px,space, electrodeTable.ChannelName(relevantElectrodeIndexMacro(CorrMacro)), format_eps);
    
    saveas(gcf,figname_fig);
    saveas(gcf,figname_jpg);
    saveas(gcf,figname_eps);
    close all
   %% Get rms amplitude
  
  %
  resp_raw_Right=mean_raw_LFPmacro_Right(2510:3010);
  resp_norm_Right=mean_norm_LFPmacro_Right((2510:3010));
  
  rms_raw_Right=rms(resp_raw_Right);
  rms_norm_Right=rms(resp_norm_Right);

  resp_raw_Left=mean_raw_LFPmacro_Left(2510:3010);
  resp_norm_Left=mean_norm_LFPmacro_Left((2510:3010));
  
  rms_raw_Left=rms(resp_raw_Left);
  rms_norm_Left=rms(resp_norm_Left);
  
  %% Rename
  
  ipsilateral_side="Right"
  contralateral_side="Left"
  
  mean_raw_LFPmacro_Ipsi=mean_raw_LFPmacro_Right;
  mean_norm_LFPmacro_Ipsi=mean_norm_LFPmacro_Right;
  mean_raw_LFPmacro_Contra=mean_raw_LFPmacro_Left;
  mean_norm_LFPmacro_Contra=mean_norm_LFPmacro_Left;
 
  resp_raw_Ipsi=resp_raw_Right;
  resp_norm_Ipsi=resp_norm_Right;
  rms_raw_Ipsi=rms_raw_Right;
  rms_norm_Ipsi=rms_norm_Right;
  resp_raw_Contra=resp_raw_Left;
  resp_norm_Contra=resp_norm_Left;
  rms_raw_Contra=rms_raw_Left;
  rms_norm_Contra=rms_norm_Left;
end
    %% SAVE THE RESULTS IN A STRUCT
    cd(ResultsFolder)
    %load(NameResulttable);
    
    ResultsTable(CorrMacro).Px=Px;%Patient
    ResultsTable(CorrMacro).Location=Location; % General location
    ResultsTable(CorrMacro).CorrMacro=CorrMacro;%Macro position on Table
    ResultsTable(CorrMacro).ElectrodeName=ElectrodeName;
    
    %General LFP
    
    ResultsTable(CorrMacro).mean_raw_LFPmacro=mean_raw_LFPmacro;
    ResultsTable(CorrMacro).mean_norm_LFPmacro=mean_norm_LFPmacro;
    ResultsTable(CorrMacro).MacroERPgood=MacroERPgood;
    
    % Signal peaks and troughs
    ResultsTable(CorrMacro).pk_1_raw=pk_1_raw;
    ResultsTable(CorrMacro).tr_raw=tr_raw;
    ResultsTable(CorrMacro).pk_2_raw=pk_2_raw;
    ResultsTable(CorrMacro).loc_pk_1_raw=loc_pk_1_raw;
    ResultsTable(CorrMacro).loc_tr_raw=loc_tr_raw;
    ResultsTable(CorrMacro).loc_pk_1_raw=loc_pk_1_raw;
    ResultsTable(CorrMacro).pk_1_norm=pk_1_norm;
    ResultsTable(CorrMacro).tr_norm=tr_norm;
    ResultsTable(CorrMacro).pk_2_norm=pk_2_norm;
    ResultsTable(CorrMacro).loc_pk_1_norm=loc_pk_1_norm;
    ResultsTable(CorrMacro).loc_tr_norm=loc_tr_norm;
    ResultsTable(CorrMacro).loc_pk_1_norm=loc_pk_1_norm;
    
    
    % Peak to trough amplitude (raw and norm)
    ResultsTable(CorrMacro).pt_raw=pt_raw;
    ResultsTable(CorrMacro).pt_norm=pt_raw;
    
    % RMS amplitude (raw and norm)
    ResultsTable(CorrMacro).rms_raw=rms_raw;
    ResultsTable(CorrMacro).rms_norm=rms_norm;
    
    %Direction
    ResultsTable(CorrMacro).isLeft=isLeft;
    ResultsTable(CorrMacro).ipsilateral_side=ipsilateral_side;
    ResultsTable(CorrMacro).contralateral_side=contralateral_side;
    
    ResultsTable(CorrMacro).mean_raw_LFPmacro_Ipsi=mean_raw_LFPmacro_Ipsi;
    ResultsTable(CorrMacro).mean_norm_LFPmacro_Ipsi=mean_norm_LFPmacro_Ipsi;
    ResultsTable(CorrMacro).mean_raw_LFPmacro_Contra=mean_raw_LFPmacro_Contra;
    ResultsTable(CorrMacro).mean_norm_LFPmacro_Contra=mean_norm_LFPmacro_Contra;
    
    ResultsTable(CorrMacro).pk_1_raw_Ipsi=pk_1_raw_Ipsi;
    ResultsTable(CorrMacro).tr_raw_Ipsi=tr_raw_Ipsi;
    ResultsTable(CorrMacro).pk_2_raw_Ipsi=pk_2_raw_Ipsi;
    ResultsTable(CorrMacro).loc_pk_1_raw_Ipsi=loc_pk_1_raw_Ipsi;
    ResultsTable(CorrMacro).loc_tr_raw_Ipsi=loc_tr_raw_Ipsi;
    ResultsTable(CorrMacro).loc_pk_2_raw_Ipsi=loc_pk_2_raw_Ipsi;
    ResultsTable(CorrMacro).pk_1_raw_Contra=pk_1_raw_Contra;
    ResultsTable(CorrMacro).tr_raw_Contra=tr_raw_Contra;
    ResultsTable(CorrMacro).pk_2_raw_Contra=pk_2_raw_Contra;
    ResultsTable(CorrMacro).loc_tr_raw_Contra=loc_tr_raw_Contra;
    ResultsTable(CorrMacro).loc_pk_2_raw_Contra=loc_pk_2_raw_Contra;
    ResultsTable(CorrMacro).pk_1_norm_Ipsi=pk_1_norm_Ipsi;
    ResultsTable(CorrMacro).tr_norm_Ipsi=tr_norm_Ipsi;
    ResultsTable(CorrMacro).pk_2_norm_Ipsi=pk_2_norm_Ipsi;
    ResultsTable(CorrMacro).loc_pk_1_norm_Ipsi=loc_pk_1_norm_Ipsi;
    ResultsTable(CorrMacro).loc_tr_norm_Ipsi=loc_tr_norm_Ipsi;
    ResultsTable(CorrMacro).loc_pk_2_norm_Ipsi=loc_pk_2_norm_Ipsi;
    ResultsTable(CorrMacro).pk_1_norm_Contra=pk_1_norm_Contra;
    ResultsTable(CorrMacro).tr_norm_Contra=tr_norm_Contra;
    ResultsTable(CorrMacro).pk_2_norm_Contra=pk_2_norm_Contra;
    ResultsTable(CorrMacro).loc_pk_1_norm_Contra=loc_pk_1_norm_Contra;
    ResultsTable(CorrMacro).loc_tr_norm_Contra=loc_tr_norm_Contra;
    ResultsTable(CorrMacro).loc_pk_2_norm_Contra=loc_pk_2_norm_Contra;
    
    ResultsTable(CorrMacro).pt_raw_Ipsi=pt_raw_Ipsi;
    ResultsTable(CorrMacro).pt_norm_Ipsi=pt_norm_Ipsi;
    ResultsTable(CorrMacro).pt_raw_Contra=pt_raw_Contra;
    ResultsTable(CorrMacro).pt_norm_Contra=pt_norm_Contra;
    ResultsTable(CorrMacro).resp_raw_Ipsi=resp_raw_Ipsi;
    ResultsTable(CorrMacro).resp_norm_Ipsi=resp_norm_Ipsi;
    ResultsTable(CorrMacro).rms_raw_Ipsi=rms_raw_Ipsi;
    ResultsTable(CorrMacro).rms_norm_Ipsi=rms_norm_Ipsi;
    ResultsTable(CorrMacro).resp_raw_Contra=resp_raw_Contra;
    ResultsTable(CorrMacro).resp_norm_Contra=resp_norm_Contra;
    ResultsTable(CorrMacro).rms_raw_Contra=rms_raw_Contra;
    ResultsTable(CorrMacro).rms_norm_Contra=rms_norm_Contra;
    
    save(NameResulttable,'ResultsTable','-v7.3')
close all
toc
end

