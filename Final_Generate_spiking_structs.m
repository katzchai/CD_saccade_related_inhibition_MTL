clear all

baseFolderPath = 'D:\Matlab\neurol\data\Saccade_New\Chaim_125_Demo_2\Saccade_ParsedData';

cd(baseFolderPath)
load(fullfile(baseFolderPath, 'electrodeTable.mat'));
load(fullfile(baseFolderPath, 'trialTable.mat'));
load(fullfile(baseFolderPath, 'trialTableImg.mat'));
rawData = load(fullfile(baseFolderPath, 'allData.mat'));

%%Patient and Information %Change the selected fields accordingly

Px={'TWH125_2'};%Change
ResultsFolder='D:\Matlab\neurol\data\Saccade_New\RsultsSaccade_v2\TWH125_2';
NameTable='Expanded_TWH125_2_Spike_epochs.mat';

%
relevantElectrodeIndexSpike = find(electrodeTable.Spike);
relevantImgOnIndices= 1:height(trialTableImg);
NumberofImgOnset=numel(relevantImgOnIndices);

%% Eliminate blinks

relevantSaccadeIndices = 1:height(trialTable);
relevantSaccadeIndices_noblinks = relevantSaccadeIndices(~isinf(trialTable.SaccadeStartStop(:,1)) & ~isnan(trialTable.SaccadeStartStop(:,1)));

NumberofSaccades=numel(relevantSaccadeIndices_noblinks);

%% Start times
saccadeStartTimes = trialTable.SaccadeStartStop(relevantSaccadeIndices_noblinks, 1);

ImgOnStartTimes = trialTableImg.ImgOn(relevantImgOnIndices);
%%
timeBeforeSaccade =2500;
timeafterSaccade = 2500;

timeBeforeImgOnset =2500;
timeafterImgOnset = 2500;
%%
%% Spike Info Control

%Create Random timestimes for control (only for plotting purposes)
minTime = 1;
maxTime = numel(rawData.time);

randomTimeIdx = round((maxTime-minTime).*rand(NumberofSaccades,1) + minTime);
randomTime=rawData.time(randomTimeIdx);
%Epoch Data
epochedDataRand=createEpochsWithDuration(rawData.data, rawData.time, randomTime, timeBeforeSaccade,...
    timeafterSaccade, relevantElectrodeIndexSpike);
totalSpikes_Control = squeeze(sum(epochedDataRand, 1));

%ImgOn
randomTimeIdxImgOn = round((maxTime-minTime).*rand(NumberofImgOnset,1) + minTime);
randomTimeImgOn=rawData.time(randomTimeIdxImgOn);
%Epoch Data
epochedDataRand_ImgOn=createEpochsWithDuration(rawData.data, rawData.time, randomTimeImgOn, timeBeforeImgOnset,...
    timeafterImgOnset, relevantElectrodeIndexSpike);
totalSpikes_Control_ImgOn = squeeze(sum(epochedDataRand_ImgOn, 1));
%%
%% Spike analysis Saccade Periods

%Epoch data saccades
epochedDataSaccade= createEpochsWithDuration(rawData.data, rawData.time, saccadeStartTimes, timeBeforeSaccade,...
    timeafterSaccade, relevantElectrodeIndexSpike);
totalSpikes_Saccade = squeeze(sum(epochedDataSaccade, 1));

%% Epoch data Image Onset and ITI onset
epochedDataImgOn= createEpochsWithDuration(rawData.data, rawData.time, ImgOnStartTimes, timeBeforeImgOnset,...
    timeafterImgOnset, relevantElectrodeIndexSpike);
totalSpikes_ImgOn = squeeze(sum(epochedDataImgOn, 1));

epochedDataITIOn= createEpochsWithDuration(rawData.data, rawData.time, ImgOnStartTimes+4000, timeBeforeImgOnset,...
    timeafterImgOnset, relevantElectrodeIndexSpike);
totalSpikes_ITIOn = squeeze(sum(epochedDataITIOn, 1));
%%
cd(ResultsFolder)
save(NameTable,'epochedDataRand','totalSpikes_Control','epochedDataRand_ImgOn','totalSpikes_Control_ImgOn',...
    'epochedDataSaccade','totalSpikes_Saccade','epochedDataImgOn','totalSpikes_ImgOn',...
    'epochedDataITIOn','totalSpikes_ITIOn','-v7.3')

