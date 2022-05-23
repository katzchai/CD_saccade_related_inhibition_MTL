% 
% % Load the data
 clear all

baseFolderPath = 'D:\Matlab\neurol\data\Saccade_New\Chaim_125_Demo_1\Saccade_ParsedData';%where the data is located
% 
cd (baseFolderPath)
load(fullfile(baseFolderPath, 'electrodeTable.mat')); %Table with macro, micro and neuron information *Note 1.1 in Word document
load(fullfile(baseFolderPath, 'trialTable.mat')); %Table with trial information: ImageNumber, saccades starts and stops and duration

rawData = load(fullfile(baseFolderPath, 'allData.mat'));%all data


%Get direction array
load('D:\Matlab\neurol\data\Saccade_New\Chaim_125_Demo_1\TWH125OneNew.mat')%provided by Chaim
directionVector=directionVectorCleanedNew;%


%Change for every patient
ResultsFolder='D:\Matlab\neurol\data\Saccade_New\RsultsSaccade_v3\TWH125_1';%where struct will be located
GeneralResultsFolder='D:\Matlab\neurol\data\Saccade_New\RsultsSaccade_v3\TWH125_1\Figures\LFP\General\';%folders for figures
DirectionResultsFolder='D:\Matlab\neurol\data\Saccade_New\RsultsSaccade_v3\TWH125_1\Figures\LFP\Direction\';%folders for figures
NameResulttable='E_TWH125_ResultsTable_New_LFP_OnlyShort500ms.mat';%Table of results
load('D:\Matlab\neurol\data\Saccade_New\RsultsSaccade_v2\TWH125_1\corresp_electrode.mat');%This is a table I made
% %to be able to get a general electrode name (e.g. LAM), macroelectrode
% %number, micro electrode number and if the electrode was relevant for
% %anlysis (e.g. corresponds to MTL). Each row is a neuron. 


 %% We are analyzing only the first macro electrode
 Px={'TWH125_1'};%Change %struct
 Location=corresp_electrode(:,1);
 macros=str2double(corresp_electrode(:,2));
 elecs_needed=unique (macros)
 
%% Patient and Information %Change the selected fields accordingly
for ii=1:size(elecs_needed)

 CorrMacro=elecs_needed(ii)
Temploc= find(macros == CorrMacro);
Temploc=Temploc(1)
electrodename= corresp_electrode(Temploc,1)
Relevant_loc= corresp_electrode(Temploc,4)
%% If location does not belong to Mesial structure--> terminate
if Relevant_loc == 'N' %only get electrodes from MTL or OCC
    Note= 'Location not Relevant '
    StructFields_ifFail_LFP  
elseif Relevant_loc == 'Y'
%%
ResultsTable=F_MakeResultsTableNew_LFPonlyLong(Px,ResultsFolder,NameResulttable,Location,rawData,...
    electrodeTable,trialTable,CorrMacro,directionVector)
end 
close all
end