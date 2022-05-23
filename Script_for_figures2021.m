%%%%%%This script intends to explain and document the code generated for the
%figures of saccade paper 2021

%Figures 2A to 2D are very simple, no stats involved so they are not
%included here

%% Note 1.1 Figure 2E

% I load 2 arrays Change_decrease_all_MTL and Change_increse_all_MTL. The different 
%arrays are expressed in percent change of the FR. For this I divided the FR in the saccade trials by the mean 
%firing rate (mean of distributions of FR of 400ms duration ‘averageFiringRates’). 
%The arrays consist then of cells x change value.

% I then get the mean and standard error of each array

numel_mean_Change_decrease_all_MTL=size(Change_decrease_all_MTL,1)
mean_Change_decrease_all_MTL=mean(Change_decrease_all_MTL)
SEM_Change_decrease_all_MTL = std(Change_decrease_all_MTL)/sqrt(numel_mean_Change_decrease_all_MTL); 

x = linspace(-2500,2500,100);
shadedErrorBar(x,mean_Change_decrease_all_MTL,SEM_Change_decrease_all_MTL,'y')
hold on

numel_mean_Change_increase_all_MTL=size(Change_increase_all_MTL,1)
mean_Change_increase_all_MTL=mean(Change_increase_all_MTL)
SEM_Change_increase_all_MTL = std(Change_increase_all_MTL)/sqrt(numel_mean_Change_increase_all_MTL); 

shadedErrorBar(x,mean_Change_decrease_all_MTL,SEM_Change_decrease_all_MTL,'r')

xlim([-500 500])

%Insert of figure 2E
% Note 1.2. Latencies
%For that I used the Probability Density ('PDF_Saccade') that consists of an array from 1 x 100 samples. This is then not
%in ms exactly. So I used that array and interpolated using imresize to
%5000 to get an array in ms. 

%%%%%%%%%%%%Need to be discussed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resized_PDF_Per_neuron_decrease=imresize(PDF_Per_neuron_decrease,[31,5000]);

minval=[]
for ii=1:size(resized_PDF_Per_neuron_decrease,1)
[l(ii),m(ii)]=min(resized_PDF_Per_neuron_decrease(ii))
end

min_vals_PDF_decrease=m-2500;% to allign to zero

%then repeated the same for increases


resized_PDF_Per_neuron_decrease=imresize(PDF_Per_neuron_decrease,[31,5000]);

minval=[]
for ii=1:size(resized_PDF_Per_neuron_decrease,1)
[l(ii),m(ii)]=min(resized_PDF_Per_neuron_decrease(ii))
end

min_vals_PDF_decrease=m-2500;% to allign to zero

%% Figure 2F. LFP analysis
% First for the LFP I ran a couple of scripts, similar to the cell by cell
% analysis, but in this case for the electrodes of interest (also included
% in the folder, "FinalCreateResultsTable_LFPonly" that calls the
% Final_MakeResultsTableNew_LFPonly.

% In general for the ERPs we took the LFP from the 1st macro electrode and
% converted to V (using the NLX file values), clean from 60Hz, and
% normalize by the std of the signal during the entire experiment. Once I
% had the LFP's from all relevant electrodes. All the ERP's from those
% electrodes in one subject were averaged. When I got the values of those
% ERP's I ploted the mean plust SE of the ERP's from all subjects.

Mean_Comp_ERPs_MTL=mean(ERP_MTL);% where ERP_MTL is a matrix of the 10 (subjects) x 5001 samples (ms)
SEM_Comp_ERPs_MTL = std(ERP_MTL)/sqrt(size(ERP_MTL,1));

f=figure('Name','ERPS','Color','w','Position',[800 100 900 900])

x=-500:1:500
hi=shadedErrorBar(x,Mean_Comp_ERPs_MTL,50,SEM_Comp_ERPs_MTL,50,'m')

line1 = lineplot(0, 'v', '-g', 'LineWidth', 1);
xlim([-500 500])
box off
title('ERP MTL')
ylabel('Norm amplitude (a.u.)')
xlabel('Time(ms)')

%% Figure 3
%For cell typing Kramyay would have to provide the code.
%this figure only consisted in detecting the cells that where modulated
%with the cell typing provided by Kramay. 
%For the waveforms I just did:


figure;
for ii=1:size (Type_dec_MTL,1)
    plot(normalize(WF_dec_MTL(:,ii)), 'color',[0 1 Type_dec_MTL(ii)],'LineWidth',2)
    hold on
end

title ('Decrease MTL')


%Where WF WF_dec_MTL are the 31 modulated cells in the MTL by 256 points of
%the WF
%And Type_dec_MTL is just an array of O(narrow) and 1 (broad) for the 31
%neurons. 

%And the same for increases

%% Figure 4. Image onset. A is an example cell. This has already been explained in the unit by unit analysis
% and the other subfigures are just counts. 

%% Figure 5. 
% Figure 5C. As in the first ERP figure (but separated by direction).
% First got the mean ERP's per subject and then got the mean ERP for all
% subjects, therefore the n=10 (per side). 

% Figure 5D. and E
% This analysis was done also by direction and by subject. So I took the
% mean CHANGE of the firing rate (mean FR during the 400 ms saccade period/
% mean FR of the Bootsrapped control periods). This was done only for the
% neurons that were modulated (Decreases only, in the contralateral or ipsilateral
% side. 
% Then I took the mean by subject and ended up with 10 values.

% The RMS amplitude was done as explained above, and then I took all the
% RMS values per subject to have an average per subject (10 ipsilateral and
% 10 contralateral).

% Then I conducted mean and SE and did stats as follows:
 mean_RMS_Amplitude_Ipsi=mean(RMS_MTL_Ipsi);
  SE_RMS_Ipsi=std(RMS_MTL_Ipsi)/sqrt(size(RMS_MTL_Ipsi,1));
 mean_RMS_Amplitude_Contra=mean(RMS_MTL_Contra);
  SE_RMS_Contra=std(RMS_MTL_Contra)/sqrt(size(RMS_MTL_Contra,1));
 
  [h,p,ci,stats] = ttest(RMS_MTL_Ipsi,RMS_MTL_Contra)
 figure;
  subplot (1,2,1)
b=bar([1:2],[mean_RMS_Amplitude_Ipsi mean_RMS_Amplitude_Contra])
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 1];b.CData(2,:) = [1 0 0];
hold on
errorbar([1:2],[mean_RMS_Amplitude_Ipsi mean_RMS_Amplitude_Contra],...
    [SE_RMS_MTL_Ipsi SE_RMS_MTL_Contra],'Color','k','LineStyle','none')
set(gca, 'XTickLabel', {'Ipsi' 'Contra'})
ylabel('Normalized RMS Amplitude (a.u.)')

ax=gca
ax.FontSize=20
box off
 

mean_FR_Saccade_DI_Ipsi=mean(FR_Saccade_DI_Ipsi);
 SE_FR_Saccade_DI_Ipsi=std(FR_Saccade_DI_Ipsi)/sqrt(size(FR_Saccade_DI_Ipsi,1));
mean_FR_Saccade_DC_Contra=mean(FR_Saccade_Per_DC_MTL_Contra);
 SE_FR_Saccade_DC_Contra=std(FR_Saccade_Per_DC_MTL_Contra)/sqrt(size(FR_Saccade_Per_DC_MTL_Contra,1));

 
 [h,p,ci,stats] = ttest(FR_Saccade_DI_Ipsi,FR_Saccade_Per_DC_MTL_Contra)

 
subplot (1,2,2)
b=bar([1:2],[mean_FR_Saccade_DI_Ipsi mean_FR_Saccade_DC_Contra])
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 1];b.CData(2,:) = [1 0 0];
hold on
errorbar([1:2],[mean_FR_Saccade_DI_Ipsi mean_FR_Saccade_DC_Contra],...
    [SE_FR_Saccade_DI_Ipsi SE_FR_Saccade_DC_Contra],'Color','k','LineStyle','none')
set(gca, 'XTickLabel', {'Ipsi' 'Contra'})
ylabel('% change in FR')

ax=gca
ax.FontSize=20
box off

% For the correlation plot:


 ALL_RMS=[RMS_MTL_Ipsi; RMS_MTL_Contra];

 
 plot (All_decreases,ALL_RMS,'.')
 
[RHO, pval]=corr(All_decreases,ALL_RMS,'type','Spearman')
lsline
title ('Change in FR/RMS amplitude')

ax=gca
ax.FontSize=20
 
 xlabel('FR percent decrease','FontSize',20)
 ylabel('ERP Amplitude','FontSize',20)

 hold on
 
 plot(FR_Saccade_DI_Ipsi,RMS_MTL_Ipsi,'bo')
 plot(FR_Saccade_Per_DC_MTL_Contra,RMS_MTL_Contra,'ro')
 
 box off