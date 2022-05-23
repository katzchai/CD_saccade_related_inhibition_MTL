clear all

baseFolderPath = 'D:\Matlab\neurol\data\Saccade_New\Chaim_125_Demo_2\Saccade_ParsedData';%Change

cd(baseFolderPath)
load(fullfile(baseFolderPath, 'electrodeTable.mat'));
load(fullfile(baseFolderPath, 'trialTable.mat'));
rawData = load(fullfile(baseFolderPath, 'allData.mat'));

%%Patient and Information %Change the selected fields accordingly


ResultsFolder='D:\Matlab\neurol\data\Saccade_New\RsultsSaccade_v2\TWH125_2';%Change
NameFRtable='TWH125_2_FR_New_IMG.mat';%Change

%
relevantElectrodeIndexSpike = find(electrodeTable.Spike);
relevantSaccadeIndices = 1:height(trialTable);
NumberofSaccades=40;

%%

time_bef_cont=750;
time_aft_cont=750;

%Create Random timestimes for control
minTime = 1;
maxTime = numel(rawData.time);
clear averageFiringRates
clear totalSpikeCounts
tic

for iteration = 1:1070
        randomTimeIdx = round((maxTime-minTime).*rand(NumberofSaccades,1) + minTime);
        %r_range = [min(randomTimeIdx) max(randomTimeIdx)]
        randomTime=rawData.time(randomTimeIdx);
        %Epoch Data
        epochedDataRand=createEpochsWithDuration(rawData.data, rawData.time, randomTime, time_bef_cont,...
            time_aft_cont, relevantElectrodeIndexSpike);
        totalSpikes = squeeze(sum(epochedDataRand, 1));
        for i = 1:size(totalSpikes, 1)
            totalSpikesBinned(i,:) = histcounts(find(totalSpikes(i,:)), 0:50:1500);
        end
        
        averageFiringRatesImg(iteration, :) = sum(totalSpikesBinned, 2)./NumberofSaccades;        
   toc     
end
cd(ResultsFolder)
save(NameFRtable,'averageFiringRatesImg')

