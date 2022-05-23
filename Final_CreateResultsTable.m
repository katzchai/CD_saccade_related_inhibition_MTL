% % % 
%%%%%%% This script prepares the data necessary to run the neuron-by-neuron
%%%%%%% analysis for saccadic modulation. This is done by subject. 
% % % Load the data
% clear all
% 
% baseFolderPath = 'D:\Matlab\neurol\data\Saccade_New\Chaim_090_Demo\Saccade_ParsedData';% folder where the data  of that specific subject is located
% 
% cd (baseFolderPath)
% load(fullfile(baseFolderPath, 'electrodeTable.mat')); %Table with macro, micro and neuron information *Note 1.1 in Word document
% load(fullfile(baseFolderPath, 'trialTable.mat')); %Table with trial information: ImageNumber, saccades starts and stops and duration
% %%%% %*Note 1.2 in Word document
% rawData = load(fullfile(baseFolderPath, 'allData.mat'));%data with all neural information according to the electrodeTable
% load(fullfile(baseFolderPath, 'trialTableImg.mat'));%Table with times of image onset only (40 trials always)
% 
% 
% %%%% %Get direction array 
% load('D:\Matlab\neurol\data\Saccade_New\Chaim_090_Demo\TWH090OneNew.mat')%load direction vector (provided by Chaim)
% directionVector=directionVectorCleanedNew;%just naming the direction vector
% 
% 
% %%%%%%% %Change for every patient
% ResultsFolder='D:\Matlab\neurol\data\Saccade_New\RsultsSaccade_v3\TWH090';%Folder to save results
% GeneralResultsFolder='D:\Matlab\neurol\data\Saccade_New\RsultsSaccade_v3\Figures\General\';%Folder for figures of general results
% DirectionResultsFolder='D:\Matlab\neurol\data\Saccade_New\RsultsSaccade_v3\Figures\Direction\';%Folder for figures of directional reuslts
% TableFR=load('D:\Matlab\neurol\data\Saccade_New\RsultsSaccade_v2\TWH090\TWH090_FR_New.mat');%Table of Firing rates of randomized periods
% %during the entire experiment. This was created by Kramay: consists of 1000 distributions of FR of 400ms duration by the number of saccades. Name:
% %%%%%%%% %Bootstrapped_sampling_FR_New *Note 1.3 in Word document.
% TableFR_IMG=load('D:\Matlab\neurol\data\Saccade_New\RsultsSaccade_v2\TWH090\TWH090_FR_New_IMG.mat');%Same as TableFR is a Table of randomized periods
% %%%%%%%%% %during the entire experiment. Consists of 1000 distributions of FR of 1500ms duration by the number of trials (40).
% load('D:\Matlab\neurol\data\Saccade_New\RsultsSaccade_v2\TWH090\corresp_electrode.mat');%This is a table I made
% %%%%%%%%% %to be able to get a general electrode name (e.g. LAM), macroelectrode
% %%%%%%%%% %number, micro electrode number and if the electrode was relevant for
% %%%%%%%%% %anlysis (e.g. corresponds to MTL). Each row is a neuron. 
% TableSpikeTrains=load('D:\Matlab\neurol\data\Saccade_New\RsultsSaccade_v2\TWH090\Expanded_TWH090_Spike_epochs.mat');%This table has
% %%%%%%%%% %epoched data of spike trains. *Note 1.4 in word document
% SpikeWF=load('D:\Matlab\neurol\data\Saccade_New\Chaim_090_Demo\Saccade_ParsedData\spikeWF.mat');%this has information of the waveform per nueron


 %%
NameResulttable="After_code_revisions"
%% Patient and Information %Change the selected fields accordingly
for ii= 1:89;
    Neuron=ii;
   disp (Neuron);
    current_corr_electrode=corresp_electrode(Neuron,:);%using the corresp_electrode table I created
Px={'TWH090'};%Change

CorrMacro=str2double(current_corr_electrode(2));% M
CorrMicro=str2double(current_corr_electrode(3));%
Location=current_corr_electrode{1};%Change: LAH/PHC/LAM/CTX/CUN... %struct
Relevant_loc=current_corr_electrode{4}; %Change: If location is not a Mesial Structure or Occipital it was not analyzed 

%% If location does not belong to Mesial structure--> terminate
if Relevant_loc == 'N'
    Note= 'Location not Relevant '
    StructFields_ifFail  
elseif Relevant_loc == 'Y'
%%
ResultsTable=FinalFinal_MakeResultsTableNew(Px,ResultsFolder,NameResulttable,Location,Relevant_loc,SpikeWF,...
    electrodeTable,trialTable,trialTableImg,Neuron,CorrMacro,CorrMicro,directionVector,TableFR,TableFR_IMG,TableSpikeTrains)
end 
close all
end