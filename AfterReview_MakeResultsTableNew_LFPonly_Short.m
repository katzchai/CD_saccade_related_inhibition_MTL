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
ElectrodeName=electrodeTable.ChannelName(relevantElectrodeIndexMacro(CorrMacro))
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

 
    %% SAVE THE RESULTS IN A STRUCT
    cd(ResultsFolder)
    %load(NameResulttable);
    
    ResultsTable(CorrMacro).Px=Px;%Patient
    ResultsTable(CorrMacro).Location=Location; % General location
    ResultsTable(CorrMacro).CorrMacro=CorrMacro;%Macro position on Table
    ResultsTable(CorrMacro).ElectrodeName=ElectrodeName;
    ResultsTable(CorrMacro).MacroERPPosNeg=MacroERPPosNeg;%1-Neg 0-Pos
    
    %General LFP
    
    ResultsTable(CorrMacro).mean_raw_LFPmacro=mean_raw_LFPmacro;
    ResultsTable(CorrMacro).mean_norm_LFPmacro=mean_norm_LFPmacro;
    ResultsTable(CorrMacro).MacroERPgood=MacroERPgood;

    
    % RMS amplitude (raw and norm)
    ResultsTable(CorrMacro).rms_raw=rms_raw;
    ResultsTable(CorrMacro).rms_norm=rms_norm;
    
   
    
    save(NameResulttable,'ResultsTable','-v7.3')
close all
toc
end

