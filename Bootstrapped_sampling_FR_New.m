clear all

baseFolderPath = 'D:\Matlab\neurol\data\Saccade_New\Chaim_089_Demo\Saccade_ParsedData';

cd(baseFolderPath)
load(fullfile(baseFolderPath, 'electrodeTable.mat'));
load(fullfile(baseFolderPath, 'trialTable.mat'));
rawData = load(fullfile(baseFolderPath, 'allData.mat'));

%%Patient and Information %Change the selected fields accordingly

Px={'TWH89'};%Change
ResultsFolder='D:\Matlab\neurol\data\Saccade_New\RsultsSaccade_v2\TWH89';
NameFRtable='TWH89_FR_New_IMG.mat';

%
relevantElectrodeIndexSpike = find(electrodeTable.Spike);
relevantSaccadeIndices = 1:height(trialTable);
NumberofSaccades=numel(relevantSaccadeIndices);

%%

time_bef_cont=200;
time_aft_cont=200;

%Create Random timestimes for control
minTime = 1;
maxTime = numel(rawData.time);
clear averageFiringRates
clear totalSpikeCounts
tic

for iteration = 1:1000
        randomTimeIdx = round((maxTime-minTime).*rand(NumberofSaccades,1) + minTime);
        %r_range = [min(randomTimeIdx) max(randomTimeIdx)]
        randomTime=rawData.time(randomTimeIdx);
        %Epoch Data
        epochedDataRand=createEpochsWithDuration(rawData.data, rawData.time, randomTime, time_bef_cont,...
            time_aft_cont, relevantElectrodeIndexSpike);
        totalSpikes = squeeze(sum(epochedDataRand, 1));
        for i = 1:size(totalSpikes, 1)
            totalSpikesBinned(i,:) = histcounts(find(totalSpikes(i,:)), 0:50:400);
        end
        
        [max_point, max_point_idx] = (max(totalSpikesBinned, [],2));
        all_max_points(iteration, :) = max_point'./NumberofSaccades;
%         all_max_points_idx(iteration, :) = max_point_idx'.*25;
        
        [min_point, min_point_idx] = (min(totalSpikesBinned, [],2));
        all_min_points(iteration, :) = min_point'./NumberofSaccades;
        
        averageFiringRates(iteration, :) = sum(totalSpikesBinned, 2)./NumberofSaccades;        
   toc     
end
cd(ResultsFolder)
save(NameFRtable,'all_max_points','all_min_points','averageFiringRates')

