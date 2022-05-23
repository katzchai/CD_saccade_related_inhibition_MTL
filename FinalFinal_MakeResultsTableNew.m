
%This function was runned for analysis per neuron to evaluate if there was
%saccadic or image onset modulation. The results were saved in a struct
%that was later used to combine information from all neurons and all
%subjects. 

function [ResultsTable]=FinalFinal_MakeResultsTableNew(Px,ResultsFolder,NameResulttable,Location,Relevant_loc,SpikeWF,...
    electrodeTable,trialTable,trialTableImg,Neuron,CorrMacro,CorrMicro,directionVector,TableFR,TableFR_IMG,TableSpikeTrains)
    
%% First see if neuron has adequate FR
% This evaluatates if the mean FR is Higher than 0.07Hz per second,
% otherwhise this neuron was discarted
tic

FR_50perc = mean(TableFR.averageFiringRates(:,Neuron))/.4
% If firing rate is too low --> terminate
if FR_50perc <= 0.07
    FRGood=0;
    Note=  'Firing Rate too low ' %Notes
%     StructFields_ifFail_v2
    return
elseif FR_50perc > 0.07
    promptFR = 'Firing Rate is Ok '
    FRGood=1
end
%% Other specs of electrode/Neuron
%this was to find the designated name of the neuron and the position it had
%in the electrode table

% Find Specific electrode positions in the electrode Table
relevantElectrodeIndexSpike = find(electrodeTable.Spike); %in the electrode table find the sections that correspond to specific neurons

%Denominations
%Neuron
CorrNeuronPos=relevantElectrodeIndexSpike(Neuron); %possition of that neuron in the electrode table. 

%% Saccade Info (get the number of saccades without blinks)

NumberofSaccades=size(TableSpikeTrains.epochedDataSaccade,1)

%note: by this point the blinks were already eliminated to create the
%TableSpikeTrains
%% Spike Info Control

%Create Random timestimes for control (only for plotting purposes)

totalSpikes_Control_Binned = histcounts(find(TableSpikeTrains.totalSpikes_Control(Neuron,:)), 0:50:5000)/NumberofSaccades;%As explained in document
%TableSpikeTrains.totalSPikes_Control are spike trains aligned to
%randomized points during the experiments. These were created for
%visualization or quality control only

totalSpikes_Control_Binned_persec=totalSpikes_Control_Binned./.05;
FR_General_ControlPer=mean(totalSpikes_Control_Binned_persec);
FR_Control_Rand=mean(totalSpikes_Control_Binned_persec(8:16));% firing rate of 400 ms previous to the "control periods" (-200 +200 in randomize periods)

%
FR_Control_Per=totalSpikes_Control_Binned_persec(46:54);%firing rate of 400 ms in aligned to control periods (-200 +200 in randomize periods)
mean_FR_Control_Per=mean(FR_Control_Per);
SE_fr_Control= std(FR_Control_Per)/sqrt(length(FR_Control_Per));


%% Asses stability of data (*Note 2.1 in Word document)
distribution_bootstrap_control=TableFR.averageFiringRates(:,Neuron)/.4; % Distributions of FR of 400ms duration 1000 iterations. Explained in first 
%section of Word document.

percUP_Con=prctile(distribution_bootstrap_control,97.5);
percDown_Con=prctile(distribution_bootstrap_control,2.5);

f=figure('Name','Stability','Color','w','Position',[800 100 900 900])
subplot (2,4,1:2)
h = histogram(distribution_bootstrap_control);
set(h,'EdgeColor','k');
h.FaceColor=[.5 .5 .5]
box off
hold on;

line1 = lineplot(percUP_Con, 'v', '--r', 'LineWidth', 1);
line2 = lineplot(percDown_Con, 'v', '--r', 'LineWidth', 1);
lineFR1=lineplot(FR_General_ControlPer, 'v', 'b', 'LineWidth', 2);%General firing rate from the entire randomize control periods
lineFR2=lineplot(FR_Control_Rand, 'v', 'g', 'LineWidth', 2);% firing rate of 400 ms previous to the "control periods" (-200 +200 in randomize periods)
lineFR3=lineplot(mean_FR_Control_Per, 'v', 'm', 'LineWidth', 2); %firing rate of 400 ms in control periods (-200 +200 in randomize periods)
leg=legend([line1 lineFR1 lineFR2 lineFR3],{'Limit',['General  ' num2str(FR_General_ControlPer)],...
    ['Stability  ' num2str(FR_Control_Rand)],['Control  ' num2str(mean_FR_Control_Per)]},'Location','northeastoutside',...
    'box','off')
box off
%

if  FR_General_ControlPer < percDown_Con || FR_Control_Rand < percDown_Con || FR_Control_Rand < percDown_Con
   title ('CAUTION!! FR stability is not good');
elseif FR_General_ControlPer > percUP_Con || FR_Control_Rand > percUP_Con || FR_Control_Rand > percUP_Con
   title ('CAUTION!! FR stability is not good')
elseif FR_General_ControlPer > percDown_Con && FR_General_ControlPer < percUP_Con ||...
       FR_Control_Rand > percDown_Con && FR_Control_Rand < percUP_Con ||...
       FR_Control_Rand > percDown_Con && FR_Control_Rand < percUP_Con
   title ('FR stability is good');
end

% Plot waveform
subplot (2,4,3)
title ('Waveform')
hold on
plot(SpikeWF.spikeInfo.mean_wf,'k')
plot(SpikeWF.spikeInfo.mean_wf(:,Neuron),'r','lineWidth',2)

% ISI
subplot (2,4,4)
allSpikeTrains_Control = squeeze(TableSpikeTrains.epochedDataRand(:,Neuron,:));
spike_occurrences=find(allSpikeTrains_Control');
spikeIntervals = spike_occurrences(2:length(spike_occurrences)) - spike_occurrences(1:length(spike_occurrences) - 1);
binSize = 50;                                            % 1 ms bins
x = [1:binSize:5000];
intervalDist = hist(spikeIntervals(spikeIntervals < 5000), x);
intervalDist = intervalDist / sum(intervalDist) / binSize; % normalize by dividing by spike number
%bar(x, intervalDist);
histfit (intervalDist,50,'exponential')

 title ('ISI');

%Plot a raster to check stability
subplot(2,1,2)

rasterplotAGPS(spike_occurrences,NumberofSaccades,5001)

box off

%% Gave a rating of confidence in that specific neuron *a manual clasification was made of the quality of the neuron
PromptGenAdj= 'General overall rating: 1 to 5 or 888 to abort    '
GenAdj=input(PromptGenAdj);%struct 

if GenAdj == 888
    Note= 'Data not good '
     StructFields_ifFail_v2 
    return
elseif GenAdj < 6
    promptGeneral = 'Good to continue '
end

    close all
%% Spike analysis Saccade Periods
%Here we calculate the FR and spike trains for the saccade periods.


totalSpikes_Saccade_Binned = histcounts(find(TableSpikeTrains.totalSpikes_Saccade(Neuron,:)), 0:50:5000)...
    /NumberofSaccades;%Spike trains in saccade trials

totalSpikes_Saccade_Binned_persec=totalSpikes_Saccade_Binned./.05;%spike trains in time

%Saccade Values *Note 2.2
FR_Saccade_Per=totalSpikes_Saccade_Binned_persec(46:54); %%%%%calculations of spike train during the saccade period%%%%%%%%%
mean_FR_Saccade_Per=mean(FR_Saccade_Per);%mean firing rate during 
SE_fr_Saccade= std(FR_Saccade_Per)/sqrt(length(FR_Saccade_Per));%standard error

%% Plot General saccade Results (*Note 2.3)

x = linspace(-2500,2500,100);%Time

halfTrialTime =2500;%Trial half time (or saccade onset)
%control
allSpikeTrains_Control = squeeze(TableSpikeTrains.epochedDataRand(:,Neuron,:));%Spike trains aligned to n(number of saccades)
%randomized points (for visualization purposes only)////Trials x Neuron x
%Time

%PDF calculations *Note 2.4
[~, colSpike_Control] =  find(allSpikeTrains_Control);% trials by time
timeStamps_ms_Control = colSpike_Control-halfTrialTime;% find timestamps where neuron spiked
fitDist_con = fitdist(timeStamps_ms_Control,'Kernel','BandWidth',50)% to do de Probability Densitity you need the fitdist first
PDF_Con = pdf(fitDist_con,x);%Probability density 

allSpikeTrains_Saccade = squeeze(TableSpikeTrains.epochedDataSaccade(:,Neuron,:));% get spike trains alligned to saccade for that 
%specific neuron
[t, colSpike_Saccade] =  find(allSpikeTrains_Saccade);% trials by time
timeStamps_ms_Saccade = colSpike_Saccade-halfTrialTime;% find timestamps where neuron spiked
fitDist_Saccade = fitdist(timeStamps_ms_Saccade,'Kernel','BandWidth',50)% to do de Probability Densitity you need the fitdist first
PDF_Saccade = pdf(fitDist_Saccade,x);%struct done

%for ploting min and max points
prob_min_Con=min(PDF_Con(46:54))
prob_min_idx_Con=find(PDF_Con == prob_min_Con);
prob_max_Con=max(PDF_Con(46:54));
prob_max_idx_Con=find(PDF_Con == prob_max_Con);

prob_min_Saccade=min(PDF_Saccade(46:54)); 
prob_min_idx_Saccade=find(PDF_Saccade == prob_min_Saccade)
prob_max_Saccade=max(PDF_Saccade(46:54));
prob_max_idx_Saccade=find(PDF_Saccade == prob_max_Saccade);

percUP=prctile(distribution_bootstrap_control,97.5); %limits of significance
percDown=prctile(distribution_bootstrap_control,2.5); %limits of significance

%figure
left_color = [0.3, 0.7, 0.9];
right_color = [0 0 0];
fig=figure('Name','GeneralFR','Color','w','Position',[1000 200 800 800])
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

subplot(4,2,1)
%Randomized periods for visualization
yyaxis left
b=bar (x,totalSpikes_Control_Binned_persec,'BarWidth', 1)
ylabel('spikes/sec')
hold on
line1 = lineplot(percUP, 'h', '--r', 'LineWidth', 1);
line2 = lineplot(percDown, 'h', '--r', 'LineWidth', 1);
ylim([0 max(totalSpikes_Control_Binned_persec)+std(totalSpikes_Control_Binned_persec)*2])
 
yyaxis right
plot (x,PDF_Con,'k','LineWidth',2)
hold on
plot (x(prob_min_idx_Con),prob_min_Con,'yo','LineWidth',2)
plot (x(prob_max_idx_Con),prob_max_Con,'ro','LineWidth',2)
xlim([-950 950])
 ylim([0 .002])
l1 = lineplot(0, 'v', 'g', 'LineWidth', 2);
l2 = lineplot(-200, 'v', '--g', 'LineWidth', 1);
l3= lineplot(200, 'v', '--g', 'LineWidth', 1);
ylabel('Probability Density')
title('Randomized Periods')

subplot(4,2,2)
%Saccade periods
yyaxis left
b=bar (x,totalSpikes_Saccade_Binned_persec,'BarWidth', 1)
ylabel('spikes/sec')
hold on
line1 = lineplot(percUP, 'h', '--r', 'LineWidth', 1);
line2 = lineplot(percDown, 'h', '--r', 'LineWidth', 1);
ylim([0 max(totalSpikes_Saccade_Binned_persec)+std(totalSpikes_Saccade_Binned_persec)*2])
 
yyaxis right
plot (x,PDF_Saccade,'k','LineWidth',2)
hold on
plot (x(prob_min_idx_Saccade),prob_min_Saccade,'yo','LineWidth',2)
plot (x(prob_max_idx_Saccade),prob_max_Saccade,'ro','LineWidth',2)
xlim([-950 950])
 ylim([0 .002])
l1 = lineplot(0, 'v', 'g', 'LineWidth', 2);
l2 = lineplot(-200, 'v', '--g', 'LineWidth', 1);
l3= lineplot(200, 'v', '--g', 'LineWidth', 1);
ylabel('Probability Density')
title('Periods aligned to saccade')

%
subplot(4,2,3)
%Raster plot control
spike_occurrences=find(allSpikeTrains_Control');

rasterplotAGPS(spike_occurrences,NumberofSaccades,5001)
currentAxis = gca;
currentAxis.XTickLabel = {'-1000','-500', '0', '500', '1000'};
xlim([2000 4000])
l1 = lineplot(3000, 'v', 'g', 'LineWidth', 2);



subplot(4,2,4)
%Raster plot control
spike_occurrences=find(allSpikeTrains_Saccade');

rasterplotAGPS(spike_occurrences,NumberofSaccades,5001)
currentAxis = gca;
currentAxis.XTickLabel = {'-1000','-500', '0', '500', '1000'};
xlim([2000 4000])
l1 = lineplot(3000, 'v', 'g', 'LineWidth', 2);



% Firing rate comparisons
%saccade period -200 to 200

subplot (4,2,5) %%%%%%This section was used for the final decision of a modulated neuron
hist(distribution_bootstrap_control)

hold on;
line1 = lineplot(percUP, 'v', '--c', 'LineWidth', 1);
line2 = lineplot(percDown, 'v', '--c', 'LineWidth', 1);
lineFR1=lineplot(mean_FR_Saccade_Per, 'v', 'g', 'LineWidth', 2);

if  mean_FR_Saccade_Per < percDown 
   title ('Significant Decrease');
elseif  mean_FR_Saccade_Per > percUP
   title ('Significant Increase');
elseif mean_FR_Saccade_Per >= percDown && mean_FR_Saccade_Per <= percUP
   title ('No effect')
end

%
subplot (4,2,6)
b=bar([1:2],[mean_FR_Control_Per mean_FR_Saccade_Per])
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 1];b.CData(2,:) = [0 1 0];
hold on
errorbar([1:2],[mean_FR_Control_Per mean_FR_Saccade_Per],[SE_fr_Control SE_fr_Saccade])
set(gca, 'XTickLabel', {'Control' 'SaccPer'})
ylim([0 max([mean_FR_Control_Per mean_FR_Saccade_Per])+.1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate if there is a rebound effect  *Note 2.5

FR_Min_Per=totalSpikes_Saccade_Binned_persec(prob_min_idx_Saccade);
mean_FR_Min_Per=mean(FR_Min_Per);
SE_fr_Min=std(FR_Min_Per)/sqrt(length(FR_Min_Per));

FR_Max_Per=totalSpikes_Saccade_Binned_persec(prob_max_idx_Saccade);
mean_FR_Max_Per=mean(FR_Max_Per);
SE_fr_Max=std(FR_Max_Per)/sqrt(length(FR_Max_Per));

distribution_bootstrap_Min=TableFR.all_min_points(:,Neuron)/.05;
mean_min=mean(distribution_bootstrap_Min);
SE_min= std(distribution_bootstrap_Min)/sqrt(length(distribution_bootstrap_Min));

percUP_Min=prctile(distribution_bootstrap_Min,97.5);
percDown_Min=prctile(distribution_bootstrap_Min,2.5);

distribution_bootstrap_max=TableFR.all_max_points(:,Neuron)/.05;
mean_max=mean(distribution_bootstrap_max);
SE_max= std(distribution_bootstrap_max)/sqrt(length(distribution_bootstrap_max));

percUP_Max=prctile(distribution_bootstrap_max,97.5);
percDown_Max=prctile(distribution_bootstrap_max,2.5);


subplot (4,2,7)%MIN

hist(distribution_bootstrap_Min)

hold on;
line1 = lineplot(percUP_Min, 'v', '--c', 'LineWidth', 1);
line2 = lineplot(percDown_Min, 'v', '--c', 'LineWidth', 1);
lineFR1=lineplot(mean_FR_Min_Per, 'v', 'y', 'LineWidth', 2);

if  mean_FR_Min_Per < percDown_Min 
   title ('Significant Decrease');
elseif  mean_FR_Min_Per > percUP_Min
   title ('Significant Increase');
elseif mean_FR_Min_Per >= percDown && mean_FR_Min_Per <= percUP_Min
   title ('No effect')
end

subplot (4,2,8)%MAX
hist(distribution_bootstrap_max)

hold on;
line1 = lineplot(percUP_Max, 'v', '--c', 'LineWidth', 1);
line2 = lineplot(percDown_Max, 'v', '--c', 'LineWidth', 1);
lineFR1=lineplot(mean_FR_Max_Per, 'v', 'r', 'LineWidth', 2);

if  mean_FR_Max_Per < percDown_Max 
   title ('Significant Decrease');
elseif  mean_FR_Max_Per > percUP_Max
   title ('Significant Increase');
elseif mean_FR_Max_Per >= percDown && mean_FR_Max_Per <= percUP_Max
   title ('No effect')
end

%%
PromptGenAdj= 'Confidence in general effect results: 1 to 5 or 888 to abort    ' 
GenAdj_RES=input(PromptGenAdj);%confidence in this particular neuron in case needed later

if GenAdj_RES == 888
    Note= 'Data not good '
     StructFields_ifFail_v2 
    return
elseif GenAdj_RES < 6
    promptGeneral = 'Good to continue '
end
%%
%Evaluate if there is modulation (although this was muanually clasified it
%was based on lines 361-267
promptSignificance = 'Is there modulation:  no=0 yes=1  ';
Significance = input(promptSignificance); 

promptTypeEffect = 'Type of effect:  NA=0 decrease=1 increase=2  rebound=3 ';
TypeEffect = input(promptTypeEffect); 

promptfigureEx= 'Is this a good example for figure? no=0 yest=1  '
figureEx=input(promptfigureEx);

 close all

%% Directionality

directionVector = directionVector(~isinf(directionVector)); %Direction vector provided by Chaim without blinks

Left_sacc=find(directionVector == -1);nLeft=numel(Left_sacc); %Indexes and n's of Left saccades
Right_sacc=find(directionVector == 1);nRight=numel(Right_sacc); %Indexes and n's of Right saccades

Neuron_epochedDataSaccade= TableSpikeTrains.epochedDataSaccade(:,Neuron,:);%Saccade spike train (reapeated from above for no reason)

epochedDataSaccade_Left= Neuron_epochedDataSaccade(Left_sacc,:,:); %epoched saccade spike trains for Left saccades
epochedDataSaccade_Right= Neuron_epochedDataSaccade(Right_sacc,:,:);%epoched saccade spike trains for Right saccades

totalSpikes_Saccade_Left = squeeze(sum(epochedDataSaccade_Left, 1))'; %spike occurrences
totalSpikes_Saccade_Right = squeeze(sum(epochedDataSaccade_Right, 1))';%spike occurrences

totalSpikes_Saccade_Binned_Left = histcounts(find(totalSpikes_Saccade_Left), 0:50:5000)...
    /nLeft;%in time

totalSpikes_Saccade_Binned_Right = histcounts(find(totalSpikes_Saccade_Right), 0:50:5000)...
    /nRight;%in time

[~, colSpike_Saccade_Left] =  find(totalSpikes_Saccade_Left); %PDF calculations (already explained)
timeStamps_ms_Saccade_Left = colSpike_Saccade_Left-halfTrialTime;
fitDist_Saccade_Left = fitdist(timeStamps_ms_Saccade_Left','Kernel','BandWidth',50);
PDF_Saccade_Left = pdf(fitDist_Saccade_Left,x);

prob_min_Saccade_Left=min(PDF_Saccade_Left(46:54));
prob_min_idx_Saccade_Left=find(PDF_Saccade_Left == prob_min_Saccade_Left);%struct done
prob_max_Saccade_Left=max(PDF_Saccade_Left(46:54));
prob_max_idx_Saccade_Left=find(PDF_Saccade_Left == prob_max_Saccade_Left);%struct done

[~, colSpike_Saccade_Right] =  find(totalSpikes_Saccade_Right);
timeStamps_ms_Saccade_Right = colSpike_Saccade_Right-halfTrialTime;
fitDist_Saccade_Right = fitdist(timeStamps_ms_Saccade_Right','Kernel','BandWidth',50);
PDF_Saccade_Right = pdf(fitDist_Saccade_Right,x);

prob_min_Saccade_Right=min(PDF_Saccade_Right(46:54));
prob_min_idx_Saccade_Right=find(PDF_Saccade_Right == prob_min_Saccade_Right);
prob_max_Saccade_Right=max(PDF_Saccade_Right(46:54));
prob_max_idx_Saccade_Right=find(PDF_Saccade_Right == prob_max_Saccade_Right);


totalSpikes_Saccade_Binned_persec_Left=totalSpikes_Saccade_Binned_Left./.05; 

%Saccade Values
FR_Saccade_Per_Left=totalSpikes_Saccade_Binned_persec_Left(46:54);% calculations of FR in saccade Periods
mean_FR_Saccade_Per_Left=mean(FR_Saccade_Per_Left);
SE_fr_Saccade_Left= std(FR_Saccade_Per_Left)/sqrt(length(FR_Saccade_Per_Left));


totalSpikes_Saccade_Binned_persec_Right=totalSpikes_Saccade_Binned_Right./.05; %Binned data

%Saccade Values
FR_Saccade_Per_Right=totalSpikes_Saccade_Binned_persec_Right(46:54);% calculations of FR in saccade Periods
mean_FR_Saccade_Per_Right=mean(FR_Saccade_Per_Right);
SE_fr_Saccade_Right= std(FR_Saccade_Per_Right)/sqrt(length(FR_Saccade_Per_Right));


%% Define side
NeuronName=electrodeTable.ChannelName(relevantElectrodeIndexSpike(Neuron)) %To define the ipsi or contralateral saccades and neurons
%Based on electrode name
Left="L"
isLeft = strfind(NeuronName,'L')

%% Plot Directional saccade Results *Note 2.6

%figure
left_color = [0.3, 0.7, 0.9];
right_color = [0 0 0];
fig=figure('Name','GeneralFR','Color','w','Position',[1000 200 800 800])
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

if isLeft == 1
    subLeft=subplot(4,2,1)
    title('Ipsi')
elseif isempty (isLeft)
    subLeft=subplot(4,2,2)
    title('Contra')
end

yyaxis left
b=bar (x,totalSpikes_Saccade_Binned_persec_Left,'BarWidth', 1)
ylabel('spikes/sec')
hold on
line1 = lineplot(percUP, 'h', '--r', 'LineWidth', 1);
line2 = lineplot(percDown, 'h', '--r', 'LineWidth', 1);
 ylim([0 max(totalSpikes_Saccade_Binned_persec_Left)+std(totalSpikes_Saccade_Binned_persec_Left)*2])
yyaxis right
plot (x,PDF_Saccade_Left,'k','LineWidth',2)
hold on
plot (x(prob_min_idx_Saccade_Left),prob_min_Saccade_Left,'yo','LineWidth',2)
plot (x(prob_max_idx_Saccade_Left),prob_max_Saccade_Left,'ro','LineWidth',2)
xlim([-950 950])
 ylim([0 .002])
l1 = lineplot(0, 'v', 'g', 'LineWidth', 2);
l2 = lineplot(-200, 'v', '--g', 'LineWidth', 1);
l3= lineplot(200, 'v', '--g', 'LineWidth', 1);
ylabel('Probability Density')


if isLeft == 1
    subLeft=subplot(4,2,2)
    title('Contra')
elseif isempty (isLeft)
    subLeft=subplot(4,2,1)
    title('Ipsi')
end

yyaxis left
b=bar (x,totalSpikes_Saccade_Binned_persec_Right,'BarWidth', 1)
ylabel('spikes/sec')
hold on
line1 = lineplot(percUP, 'h', '--r', 'LineWidth', 1);
line2 = lineplot(percDown, 'h', '--r', 'LineWidth', 1);
 ylim([0 max(totalSpikes_Saccade_Binned_persec_Right)+std(totalSpikes_Saccade_Binned_persec_Right)*2])
yyaxis right
plot (x,PDF_Saccade_Right,'k','LineWidth',2)
hold on
plot (x(prob_min_idx_Saccade_Right),prob_min_Saccade_Right,'yo','LineWidth',2)
plot (x(prob_max_idx_Saccade_Right),prob_max_Saccade_Right,'ro','LineWidth',2)
xlim([-950 950])
 ylim([0 .0015])
l1 = lineplot(0, 'v', 'g', 'LineWidth', 2);
l2 = lineplot(-200, 'v', '--g', 'LineWidth', 1);
l3= lineplot(200, 'v', '--g', 'LineWidth', 1);
ylabel('Probability Density')


% Firing rate comparisons
%saccade period -200 to 200

if isLeft == 1
    subLeft=subplot(4,2,3)
elseif isempty (isLeft)
    subLeft=subplot(4,2,4)
end

hist(distribution_bootstrap_control)

hold on;
line1 = lineplot(percUP, 'v', '--c', 'LineWidth', 1);
line2 = lineplot(percDown, 'v', '--c', 'LineWidth', 1);
lineFR1=lineplot(mean_FR_Saccade_Per_Left, 'v', 'g', 'LineWidth', 2);

if  mean_FR_Saccade_Per_Left < percDown 
   title ('Significant Decrease');
elseif  mean_FR_Saccade_Per_Left > percUP
   title ('Significant Increase');
elseif mean_FR_Saccade_Per_Left >= percDown && mean_FR_Saccade_Per_Left <= percUP
   title ('No effect')
end

if isLeft == 1
    subLeft=subplot(4,2,4)
elseif isempty (isLeft)
    subLeft=subplot(4,2,3)
end

hist(distribution_bootstrap_control)

hold on;
line1 = lineplot(percUP, 'v', '--c', 'LineWidth', 1);
line2 = lineplot(percDown, 'v', '--c', 'LineWidth', 1);
lineFR1=lineplot(mean_FR_Saccade_Per_Right, 'v', 'g', 'LineWidth', 2);

if  mean_FR_Saccade_Per_Right < percDown %%%%%%%%%%%%%% Decision of directional effect was based in lines 566 to 572
   title ('Significant Decrease');
elseif  mean_FR_Saccade_Per_Right > percUP
   title ('Significant Increase');
elseif mean_FR_Saccade_Per_Right >= percDown && mean_FR_Saccade_Per_Right <= percUP
   title ('No effect')
end


%%
PromptGenAdj_Direction= 'Grade data confidence (1 to 5 for Direction results) or 888 to abort   '
GenAdj_Direction=input(PromptGenAdj_Direction);

if GenAdj_Direction == 888
    Note= 'Data not good '
    StructFields_ifFail_v2 
    return
elseif GenAdj_Direction < 6
    promptGeneral = 'Good to continue '
end


%Evaluate if there is modulation
promptSignificance_Direction = 'Is there a direction preference?  no=0 yes=1  No modulation (NA)= 88  ';
Significance_Direction = input(promptSignificance_Direction); 


promptTypeEffect_Direction = 'Type of directional effect? II = 1   IC = 2   DI= 3  DC= 4   II/DC= 5  DI/IC =6  EQUAL= 7  NoEffect (NA) = 88   ';
TypeEffect_Direction = input(promptTypeEffect_Direction); 
% II=Increase ipsi; IC=Increase contra; DI=Decrease Ipsi; DC=Decrease
% Contra; Equal= both sides are modulated (no cases found); No effect= No directional
% modulation


promptfigureEx_Direction= 'Is this a good example for figure? no=0 yest=1 '
figureEx_Direction=input(promptfigureEx_Direction);

%%
close all


%% Image onset calculations %to see if there is a visual information effect
 %% Spike analysis Image onset Periods (*Note 2.7)

% Info

relevantImgOnIndices = 1:height(trialTableImg);% Number of trials image onset

NumberofTrials=numel(relevantImgOnIndices);%total trials

timeBeforeImgOn =2500;

%%
%Create Random timestimes for control (only for plotting purposes)

totalSpikes_ControlImgOn_Binned = histcounts(find(TableSpikeTrains.totalSpikes_Control_ImgOn(Neuron,:)), 0:50:5000)/NumberofTrials;
%Similar to saccade analysis for Image onset analysis I also calculated
%control trials for visualization. Not used in this version of the script
totalSpikes_ControlImgOn_Binned_persec=totalSpikes_ControlImgOn_Binned./.05; %trials were binned by 50 ms
FR_General_ControlImgOn_Per=mean(totalSpikes_ControlImgOn_Binned_persec); %mean (not used for analysis)
FR_ControlImgOn_Rand=mean(totalSpikes_ControlImgOn_Binned_persec(8:28));%mean in a previous time (not used for anlysis)

%
FR_ControlImgOn_Per=totalSpikes_ControlImgOn_Binned_persec(55:85);%control period spike trains (not used for analysis) in same 
%time period of analysis as image onset
mean_FR_ControlImgOn_Per=mean(FR_ControlImgOn_Per);%mean control period spike trains control (not used for analysis)
SE_fr_ControlImgOn= std(FR_ControlImgOn_Per)/sqrt(length(FR_ControlImgOn_Per));%SE

%%

distribution_bootstrap_controlIMG=TableFR_IMG.averageFiringRatesImg(:,Neuron)/1.5;%The distribution in this case consists of 1000 samples of FR at
%40 random points (equivalent to 40 img onsets) and 1500ms (equal the amount of time used for the analysis of image onset modulation 200-1700ms) 

totalSpikes_ImgOn_Binned = histcounts(find(TableSpikeTrains.totalSpikes_ImgOn(Neuron,:)), 0:50:5000)...
    /NumberofTrials; %binned spike trains Image Onset

totalSpikes_ImgOn_Binned_persec=totalSpikes_ImgOn_Binned./.05;%In time
FR_General_ImgOnPer=mean(totalSpikes_ImgOn_Binned_persec);%General Firing Rate (not used for analysis only to assess stability)
FR_ImgOn_prev=mean(totalSpikes_ImgOn_Binned_persec(8:28));%Firing rate previous to image onset (not used for analysis only to assess stability)

%ImgOn Values
FR_ImgOn_Per=totalSpikes_ImgOn_Binned_persec(55:84);%specific period of analysis similar to Ruthishauser lab
mean_FR_ImgOn_Per=mean(FR_ImgOn_Per);%mean FR in that specific period %%%% this was compared to randomized control trials
% and used for analysis
SE_fr_ImgOn= std(FR_ImgOn_Per)/sqrt(length(FR_ImgOn_Per));%SE

%% PDF (not used for analysis in the case of Image Onset)

allSpikeTrains_Control_ImgOn = squeeze(TableSpikeTrains.epochedDataRand_ImgOn(:,Neuron,:));
[~, colSpike_Control_ImgOn] =  find(allSpikeTrains_Control_ImgOn);
timeStamps_ms_Control_ImgOn = colSpike_Control_ImgOn-timeBeforeImgOn;
fitDist_con_ImgOn = fitdist(timeStamps_ms_Control_ImgOn,'Kernel','BandWidth',50)
PDF_Con_ImgOn = pdf(fitDist_con_ImgOn,x);%struct done

allSpikeTrains_ImgOn = squeeze(TableSpikeTrains.epochedDataImgOn(:,Neuron,:));
[t, colSpike_ImgOn] =  find(allSpikeTrains_ImgOn);
timeStamps_ms_ImgOn = colSpike_ImgOn-timeBeforeImgOn;
fitDist_ImgOn = fitdist(timeStamps_ms_ImgOn,'Kernel','BandWidth',50)
PDF_ImgOn = pdf(fitDist_ImgOn,x);%struct done

percUP_ImgOn=prctile(distribution_bootstrap_controlIMG,97.5);%significant values
percDown_ImgOn=prctile(distribution_bootstrap_controlIMG,2.5);%significant values
%% plot image onset (Note 2.8 example figure)
x = [1:200:5000]
Array_example_Img=find(TableSpikeTrains.totalSpikes_ImgOn(Neuron,:));
% Binning data
fig=figure('Name','GeneralFR','Color','w','Position',[1000 50 900 900]);

edges_a=[0:200:5001];%different edges to get nicer plots only
edges_b=[50:200:5051];
edges_c=[-50:200:4951];

[N_a, ~] = histcounts(Array_example_Img, edges_a);
[N_b, ~] = histcounts(Array_example_Img, edges_b);
[N_c, ~] = histcounts(Array_example_Img, edges_c);

N_all_ImgOn=mean([N_a/NumberofTrials;N_b/NumberofTrials;N_c/NumberofTrials]);
N_all_FR_ImgON=N_all_ImgOn./.2;

x = linspace(-2500,2500,25);
%  % Histograms

subplot(3,3,[1 1.5])
H_sacc=bar(x,N_all_FR_ImgON,1,'FaceColor',[0 0 0],'EdgeColor',[1 1 1])

l1 = lineplot(0, 'v', 'c', 'LineWidth', 4);
line1 = lineplot(FR_50perc, 'h', 'm', 'LineWidth', 2);
line2 = lineplot(percUP_ImgOn, 'h', '--m', 'LineWidth', 1);
line3 = lineplot(percDown_ImgOn, 'h', '--m', 'LineWidth', 1);
xlim([-1000 2500])
set(gca,'fontsize', 20);
ylabel('Firing Rate (Hz)')
box off


N_all_FR_percent_ImgOn=(N_all_FR_ImgON-FR_50perc)./abs(FR_50perc).*100;%percent change  per neuron. These values were used later
%for the MUA latency modulation analysis. 
se_perc_ImgOn=(percUP_ImgOn-FR_50perc)./abs(FR_50perc).*100;

subplot (3,3,[4 4.5])
plot(x,N_all_FR_percent_ImgOn,'Linewidth',2,'Color',[0 0 0])
xlim([-1000 2500])
ylim([-850 850])
l1 = lineplot(0, 'v', 'c', 'LineWidth', 4);
line1 = lineplot(0, 'h', 'm', 'LineWidth', 2);
line1 = lineplot(se_perc_ImgOn, 'h', '--m', 'LineWidth', 1);
line1 = lineplot(-se_perc_ImgOn, 'h', '--m', 'LineWidth', 1);
box off
set(gca,'fontsize', 20);
ylabel('%Change')

subplot (3,3,[7 7.5])
spike_ocurrences_ImgOn=find(allSpikeTrains_ImgOn')
rasterplot_real(spike_ocurrences_ImgOn,NumberofTrials,5001)
currentAxis = gca;
currentAxis.XTickLabel = { '-500', '500' '1500' '2500'};
xlim([1500 5000])
l1 = lineplot(2500, 'v', 'c', 'LineWidth', 4);
set(gca,'fontsize', 20);


subplot (8,5,3.5)
title ('Waveform')
hold on
plot(SpikeWF.spikeInfo.mean_wf(:,Neuron),'k','lineWidth',2)
axis off

subplot(6,6,10)
hist(distribution_bootstrap_controlIMG)

hold on;
line1 = lineplot(percUP_ImgOn, 'v', '--k', 'LineWidth', 1);
line2 = lineplot(percDown_ImgOn, 'v', '--k', 'LineWidth', 1);
lineFR1=lineplot(mean_FR_ImgOn_Per, 'v', 'c', 'LineWidth', 2);

set(gca,'fontsize', 20);
xlim([min(distribution_bootstrap_controlIMG)-.3  max(distribution_bootstrap_controlIMG)+.3])

%%
%Evaluate if there is modulation
PromptGenAdj_IMG= 'Grade data confidence IMG data (1 to 5 for Direction results) or 888 to abort   '
GenAdj_IMG=input(PromptGenAdj_IMG);% score of confidence (only for my own reference)

if GenAdj_IMG == 888
    Note= 'Data not good '
     StructFields_ifFail_v2 
    return
elseif GenAdj_IMG < 6
    promptGeneral = 'Good to continue '
end


promptSignificance_ImgOn = 'Is there visual modulation no=0 yes=1 data too variable or bad =88  ';
Significance_ImgOn = input(promptSignificance_ImgOn); %decision of modulation based mean FR compared to distribution.

promptTypeEffect_ImgOn = 'Increase or decrease NA=0 decrease=1 increase=2 ';
TypeEffect_ImgOn = input(promptTypeEffect_ImgOn); %struct done

promptLatency_ImgOn = 'Latency of effect NA=0 early=1 late=2 ';
LatencyEffect_ImgOn = input(promptLatency_ImgOn); %struct done

promptfigureEx_ImOn= 'Is this a good example for figure? no=0 yest=1 '
figureEx_ImgOn=input(promptfigureEx_ImOn);%struct done
%%

%%

OverallQuality=sum([GenAdj GenAdj_RES GenAdj_Direction GenAdj_IMG]);%Quality score (not used for analysis only for my own reference)

OverallQuality=sum([GenAdj GenAdj_RES GenAdj_Direction GenAdj_IMG]);%struct done

FinalClassSaccadePrompt= 'NO EFFECT = 0  EFFECT = 1  '
Final_saccade_effect=input(FinalClassSaccadePrompt);%struct done

FinalClassDirectionPrompt= 'NO EFFECT = 0  EFFECT = 1  '
Final_Directon_effect=input(FinalClassDirectionPrompt);%struct done

FinalClassImgOnPrompt = 'NO EFFECT = 0  EFFECT = 1  '
Final_ImgOn_effect=input(FinalClassImgOnPrompt);%struct done

FinalClassImgOnPrompt = 'NOT Included = 0  Included = 1  '
Final_ImgOn_effect=input(FinalClassImgOnPrompt);%struct done
toc

%% Make struct for analysis %Results where saved
 cd(ResultsFolder)
    load(NameResulttable);
    ResultsTable(Neuron).Px=Px;%Patient
    ResultsTable(Neuron).Neuron=Neuron;%Neuron number in Table
    ResultsTable(Neuron).Relevant_loc=Relevant_loc; %Is the location relevant in last analysis
    ResultsTable(Neuron).Location=Location; % General location
    ResultsTable(Neuron).CorrMacro=CorrMacro;%Macro position on Table
    
    %If terminated
    ResultsTable(Neuron).Terminated=0; %Did not complete processing: Data not good or not relevant
    ResultsTable(Neuron).Note= nan; %Notes of data
    
   % General firing rate
    ResultsTable(Neuron).FRGood=FRGood; %FR above 0.1
    ResultsTable(Neuron).FR_50perc=FR_50perc;
    
    % Other specifics of this electrode/Neuron
    ResultsTable(Neuron).NeuronDen=Location;%Name of neuron
    ResultsTable(Neuron).CorrNeuronPos=CorrNeuronPos; %Neuron position on Table
    ResultsTable(Neuron).CorrMicro=CorrMicro;
    ResultsTable(Neuron).CorrMacro=CorrMacro;
    
    % Spike analysis
    ResultsTable(Neuron).NumberofSaccades=NumberofSaccades;
    ResultsTable(Neuron).totalSpikes_Control_Binned_persec=totalSpikes_Control_Binned_persec;
    ResultsTable(Neuron).FR_Control_Per=FR_Control_Per;
    ResultsTable(Neuron).mean_FR_Control_Per=mean_FR_Control_Per;
    ResultsTable(Neuron).SE_fr_Control=SE_fr_Control;
    ResultsTable(Neuron).distribution_bootstrap_control=distribution_bootstrap_control;
    ResultsTable(Neuron).percUP=percUP;%97.5 percentile
    ResultsTable(Neuron).percDown=percDown;%2.5 percentile
    ResultsTable(Neuron).totalSpikes_Saccade_Binned_persec=totalSpikes_Saccade_Binned_persec;
    ResultsTable(Neuron).FR_Saccade_Per=FR_Saccade_Per;
    ResultsTable(Neuron).mean_FR_Saccade_Per=mean_FR_Saccade_Per;
    ResultsTable(Neuron).SE_fr_Saccade=SE_fr_Saccade;
    ResultsTable(Neuron).PDF_Con=PDF_Con;
    ResultsTable(Neuron).PDF_Saccade=PDF_Saccade;
   
    
    %Rebound
    ResultsTable(Neuron).prob_min_Con=prob_min_Con;
    ResultsTable(Neuron).prob_min_idx_Con=prob_min_idx_Con;
    ResultsTable(Neuron).prob_max_Con=prob_max_Con;
    ResultsTable(Neuron).prob_max_idx_Con=prob_max_idx_Con;
    ResultsTable(Neuron).prob_min_Saccade=prob_min_Saccade;
    ResultsTable(Neuron).prob_min_idx_Saccade=prob_min_idx_Saccade;
    ResultsTable(Neuron).prob_max_Saccade=prob_max_Saccade;
    ResultsTable(Neuron).prob_max_idx_Saccade=prob_max_idx_Saccade;
    ResultsTable(Neuron).FR_Min_Per=FR_Min_Per;
    ResultsTable(Neuron).mean_FR_Min_Per=mean_FR_Min_Per;
    ResultsTable(Neuron).SE_fr_Min=SE_fr_Min;
    ResultsTable(Neuron).FR_Max_Per=FR_Max_Per;
    ResultsTable(Neuron).mean_FR_Max_Per=mean_FR_Max_Per;
    ResultsTable(Neuron).SE_fr_Max=SE_fr_Max;
    
    % Notes of general saccade analysis (first figure)
    ResultsTable(Neuron).GenAdj=GenAdj; % General initial analysis from 1 to 5 or if aborted 888
    ResultsTable(Neuron).GenAdj_RES=GenAdj_RES; % General Saccade analysis from 1 to 5 or if aborted 888
    ResultsTable(Neuron).Significance=Significance;%'Is there modulation:  no=0 yes=1 tendency=2  ';
    ResultsTable(Neuron).TypeEffect=TypeEffect; %'Type of effect:  NA=0 decrease=1 increase=2  rebound=3 ';
    ResultsTable(Neuron).figureEx=figureEx; %'Is this a good example for figure? no=0 yest=1  '
    
    
    % Directionality
    ResultsTable(Neuron).directionVector=directionVector;
    ResultsTable(Neuron).Left_sacc=Left_sacc;%Indexes of left sacc
    ResultsTable(Neuron).Right_sacc=Right_sacc;%Indexes of right sacc
    ResultsTable(Neuron).nLeft=nLeft;%#sacc
    ResultsTable(Neuron).nRight=nRight;%#sacc
    ResultsTable(Neuron).totalSpikes_Saccade_Binned_Left=totalSpikes_Saccade_Binned_Left;
    ResultsTable(Neuron).totalSpikes_Saccade_Binned_Right=totalSpikes_Saccade_Binned_Right;
    ResultsTable(Neuron).timeStamps_ms_Saccade_Left=timeStamps_ms_Saccade_Left;
    ResultsTable(Neuron).PDF_Saccade_Left=PDF_Saccade_Left;
    ResultsTable(Neuron).prob_min_Saccade_Left=prob_min_Saccade_Left;
    ResultsTable(Neuron).prob_min_idx_Saccade_Left=prob_min_idx_Saccade_Left;
    ResultsTable(Neuron).prob_max_Saccade_Left=prob_max_Saccade_Left;
    ResultsTable(Neuron).prob_max_idx_Saccade_Left=prob_max_idx_Saccade_Left;
    ResultsTable(Neuron).timeStamps_ms_Saccade_Right=timeStamps_ms_Saccade_Right;
    ResultsTable(Neuron).PDF_Saccade_Right=PDF_Saccade_Right;
    ResultsTable(Neuron).prob_min_Saccade_Right=prob_min_Saccade_Right;
    ResultsTable(Neuron).prob_min_idx_Saccade_Right=prob_min_idx_Saccade_Right;
    ResultsTable(Neuron).prob_max_Saccade_Right=prob_max_Saccade_Right;
    ResultsTable(Neuron).prob_max_idx_Saccade_Right=prob_max_idx_Saccade_Right;
    ResultsTable(Neuron).totalSpikes_Saccade_Binned_persec_Left=totalSpikes_Saccade_Binned_persec_Left;
    ResultsTable(Neuron).FR_Saccade_Per_Left=FR_Saccade_Per_Left;
    ResultsTable(Neuron).mean_FR_Saccade_Per_Left=mean_FR_Saccade_Per_Left;
    ResultsTable(Neuron).SE_fr_Saccade_Left=SE_fr_Saccade_Left;
    ResultsTable(Neuron).totalSpikes_Saccade_Binned_persec_Right=totalSpikes_Saccade_Binned_persec_Right;
    ResultsTable(Neuron).FR_Saccade_Per_Right=FR_Saccade_Per_Right;
    ResultsTable(Neuron).mean_FR_Saccade_Per_Right=mean_FR_Saccade_Per_Right;
    ResultsTable(Neuron).SE_fr_Saccade_Right=SE_fr_Saccade_Right;
    ResultsTable(Neuron).isLeft=isLeft;
    
    % Notes about analysis of direction
    ResultsTable(Neuron).GenAdj_Direction=GenAdj_Direction;% 'So far data looks good?  yes=1 needs adjustments=2 data not good= 888 (terminate)  '
    ResultsTable(Neuron).Significance_Direction=Significance_Direction;%'Is there modulation no=0 yes=1 tendency=2  '
    ResultsTable(Neuron).TypeEffect_Direction=TypeEffect_Direction;%'II = 1   IC = 2   DI= 3  DC= 4   II/DC= 5  DI/IC =6  EQUAL= 7  NoEffect (NA) = 88 ';
    ResultsTable(Neuron).figureEx_Direction=figureEx_Direction;%'Is this a good example for figure? no=0 yest=1 '
   
    % Image data
    ResultsTable(Neuron).NumberofTrials=NumberofTrials;
    ResultsTable(Neuron).totalSpikes_ControlImgOn_Binned=totalSpikes_ControlImgOn_Binned;
    ResultsTable(Neuron).FR_ControlImgOn_Per=FR_ControlImgOn_Per;
    ResultsTable(Neuron).mean_FR_ControlImgOn_Per=mean_FR_ControlImgOn_Per;
    ResultsTable(Neuron).SE_fr_ControlImgOn=SE_fr_ControlImgOn;
    ResultsTable(Neuron).distribution_bootstrap_control=distribution_bootstrap_control;
    ResultsTable(Neuron).totalSpikes_ImgOn_Binned=totalSpikes_ImgOn_Binned;
    ResultsTable(Neuron).FR_ImgOn_Per=FR_ImgOn_Per;
    ResultsTable(Neuron).mean_FR_ImgOn_Per=mean_FR_ImgOn_Per;
    ResultsTable(Neuron).SE_fr_ImgOn=SE_fr_ImgOn;
    ResultsTable(Neuron).PDF_Con_ImgOn=PDF_Con_ImgOn;
    ResultsTable(Neuron).PDF_ImgOn=PDF_ImgOn;
    
    %Image data different binning
    ResultsTable(Neuron).Array_example_Img=Array_example_Img;
    ResultsTable(Neuron).N_all_FR_ImgOn=N_all_FR_ImgON;
    ResultsTable(Neuron).N_all_FR_percent_ImgOn=N_all_FR_percent_ImgOn;
    ResultsTable(Neuron).se_perc_ImgOn=se_perc_ImgOn;
    
    %Notes of analysis ImgOn
    ResultsTable(Neuron).Significance_ImgOn=Significance_ImgOn;
    ResultsTable(Neuron).TypeEffect_ImgOn=TypeEffect_ImgOn;
    ResultsTable(Neuron).LatencyEffect_ImgOn=LatencyEffect_ImgOn;
    ResultsTable(Neuron).figureEx_ImgOn=figureEx_ImgOn;
    
    %Overall quality of results
    ResultsTable(Neuron).OverallQuality=OverallQuality;% Max of 20
    ResultsTable(Neuron).Final_saccade_effect=Final_saccade_effect;
    ResultsTable(Neuron).Final_Directon_effect=Final_Directon_effect;
    ResultsTable(Neuron).Final_ImgOn_effect=Final_ImgOn_effect;
    ResultsTable(Neuron).Included=Included;

    %
    save(NameResulttable,'ResultsTable','-v7.3')

end
