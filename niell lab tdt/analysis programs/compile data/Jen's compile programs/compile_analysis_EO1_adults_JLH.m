%function compile developmental data
clear all
close all
dbstop if error
% dbclear all % this will exit out of db mode.
[fname pname] =uiputfile('*.ps'); psfilename=fullfile(pname,fname);  %%% get ps filename
%psfilename = 'c:/test.ps';   %%% default location
if exist(psfilename,'file')==2;delete(psfilename);end %%% check for previous file


for dataset = 1:2
    
    cells=0;
    if dataset ==1
        afiles = {'D:\Jen_ephys_data\developmental_periods\Good recordings\Adults\1month_2month old\4_22_13\analysis_adult.mat',...
            'D:\Jen_ephys_data\developmental_periods\Good recordings\Adults\5_25_13_Adultmale3mo\Adult B\analysis_5_25_13_B_rec1_adult.mat',...
            'D:\Jen_ephys_data\developmental_periods\Good recordings\Adults\1month_2month old\7_18_13_adult\analysis_07_18_13_cluster_rec1.mat',...
            'D:\Jen_ephys_data\developmental_periods\Good recordings\Adults\7_19_13_adult\analysis_cluster_adult_rec2_7_19_13.mat',...
            'D:\Jen_ephys_data\developmental_periods\Good recordings\Adults\7_25_13_deep\MouseB\Rec1\analysis_7_25_13_mouseB_rec1.mat',...
            'D:\Jen_ephys_data\developmental_periods\Good recordings\Adults\7_25_13_deep\MouseB\Rec2\analysis_cluster_data_07_25_13_mouseB_adult_rec2.mat',...
            'D:\Jen_ephys_data\developmental_periods\Good recordings\Adults\1month_2month old\9_25_13\rec1\analysis_9_25_13_P31_rec1.mat',...
            'D:\Jen_ephys_data\developmental_periods\Good recordings\Adults\9_27_13\rec1\analysis_9_27_13_adult_rec1.mat',...
            'D:\Jen_ephys_data\developmental_periods\Good recordings\Adults\9_27_13\rec2\analysis_9_27_13_adult_rec2.mat',...
            'D:\Jen_ephys_data\developmental_periods\Good recordings\Adults\9_12_13\rec1\analysis_9_12_13_adult_rec1.mat',...
            'D:\Jen_ephys_data\developmental_periods\Good recordings\Adults\9_12_13\rec2\analysis_9_12_13_adult_rec2.mat',...
            'D:\Jen_ephys_data\developmental_periods\Good recordings\Adults\11_11_13\rec1\analysis_11_11_13_adult_rec1.mat',...
            'D:\Jen_ephys_data\developmental_periods\Good recordings\Adults\11_11_13\rec2\analysis_11_11_13_adult_rec2.mat',...
            'D:\Jen_ephys_data\developmental_periods\Good recordings\Adults\11_13_13\rec1\analysis_adult_11_13_13_rec1.mat',...
            'D:\Jen_ephys_data\developmental_periods\Good recordings\Adults\11_13_13\rec2\analysis_11_13_13_rec2.mat',...
            'D:\Jen_ephys_data\developmental_periods\Good recordings\Adults\11_15_13\rec1\analysis_11_15_13_adult_rec1.mat',...
            'D:\Jen_ephys_data\developmental_periods\Good recordings\Adults\11_15_13\rec2\analysis_11_15_13_adult_rec2.mat'};
    elseif dataset ==2
        afiles = {'D:\Jen_ephys_data\developmental_periods\Good recordings\EO1_EO2\4_29_13_EO1\mouseC\analysis_4_29_13_C.mat',...
            'D:\Jen_ephys_data\developmental_periods\Good recordings\EO1_EO2\8_7_13_EO1\rec1_full_clustering\analysis_8_7_13_EO1_rec1_more_strigent.mat',...
            'D:\Jen_ephys_data\developmental_periods\Good recordings\EO1_EO2\5_22_13_EO1\analysis_rec1_A_5_22_13_strict_selection.mat',...
            'D:\Jen_ephys_data\developmental_periods\Good recordings\EO1_EO2\7_17_13_EO1\Analysis_7_17_13_cluster_7_17_13_EO1.mat',...
           'D:\Jen_ephys_data\developmental_periods\Good recordings\EO1_EO2\7_17_13_EO1\Rec2\Analysis_7_17_13_rec2_EO1.mat',...
           'D:\Jen_ephys_data\developmental_periods\Good recordings\EO1_EO2\8_8_13_EO2\rec1_full_clustering\analysis.mat',...
          'D:\Jen_ephys_data\developmental_periods\Good recordings\EO1_EO2\9_9_13_EO1\rec1\analysis_9_9_13_EO1_rec1.mat',...
          'D:\Jen_ephys_data\developmental_periods\Good recordings\EO1_EO2\9_30_13_EO1\rec1\analysis_9_30_13_EO1_rec1.mat',...
          'D:\Jen_ephys_data\developmental_periods\Good recordings\EO1_EO2\9_30_13_EO1\rec2\analysis_9_30_13_EO1_rec2.mat',...
          'D:\Jen_ephys_data\developmental_periods\Good recordings\EO1_EO2\9_9_13_EO1\rec2\analysis_9_9_13_EO1_rec2.mat',...
          'D:\Jen_ephys_data\developmental_periods\Good recordings\EO1_EO2\10_1_13_EO2\rec1\analysis_10_1_13_EO2_rec1.mat'}; %%% tg
    end
    
    N =0;
    for i = 1:length(afiles)
        
        i;
        load(afiles{i});
        n_units = length(L_ratio);
        cellrange = N+1:N+n_units;
        N=N+n_units;
        
        
        
        alldata( cellrange,1:2) = cells;
        alldata( cellrange,3) = L_ratio;
        
        %%% waveform
        alldata( cellrange,4) = trough_width;
        alldata( cellrange,5) = trough2peak;
        alldata( cellrange,6) = -trough_depth./peak_height;
        alldata( cellrange,7:25)= wv';
        
 
   
        %     for c = 1:n_units;
        %         ch = alldata(c,1);
        %         cl = alldata(c,2);
        %         min_t = min(mean_wvform (:,ch : ch+3,cl),[],1);
        %         [trough_allc trig_chan] = min(min_t);
        %         trig_chan = ch+trig_chan-1;
        %         alldata(c,7) = mean_wvform(size(mean_wvform,1),trig_chan,cl)/peak_height(c);
        %
        %         t1= squeeze(event_times_all(ch,find(idx_all(ch,:) == cl)));
        %         dt = diff(t1);
        %         dt =dt(dt<.02);
        %         n=hist(dt,.001:0.002:.02);
        %         [y alldata(c,8)] = max(n);
        %         n=hist(dt,.0005:0.001:.02);
        %         alldata(c,9) = max(n(3:8))./mean(n(15:20))
        %     end;
        
        
        %         A1(cellrange,:)=bars_A1;
        %         A2(cellrange,:)=bars_A2;
        %         w(cellrange,:)=bars_w;
        %         theta(cellrange,:)=bars_theta;
        %         bspont(cellrange,:)=barspont;
        %
        %         size(rf_width)
        %
        %         rfw(cellrange,:) = rf_width*30;
        
        
        size(drift);
        
        driftA1(cellrange,:)= field2array(drift,'A1');
        driftA2(cellrange,:)=field2array(drift,'A2');
        driftB(cellrange,:)= field2array(drift,'B');
        drift_theta_w(cellrange,:)=field2array(drift,'thetawidth');
        driftspont(cellrange,:) = field2array(drift,'spont');
        
        
        %generate tuning curve in polar plot style
%         driftthetatuning(cellrange,:) = field2array(drift,'thetatuning');
%       
%         figure
%         polar(driftthetatuning{4,:});
%         hold on
%         polar (driftthetatuning{3,:});
% %         
       
        
        
        driftwpref(cellrange,:) = field2array(drift,'wpref');
        driftwbw(cellrange,:) = field2array(drift,'bw') ;
        
        driftF1F0(cellrange,:) = field2array(drift,'F1')./field2array(drift,'F0');
     
%       driftorientfreq_all(cellrange,:)=field2array(drift, 'orientfreq_all');
        
        driftlayer(cellrange,:) = field2array(drift,'layer');
        driftdsi(cellrange,:) = field2array(drift,'dsi');
       
        %size(wvform)
        wvform(cellrange,:) = wv';
        
    end
    
    %replace all NaN in F1F0 with 0
    
    ind = find(isnan(driftF1F0));
    driftF1F0(ind)=0;
   
    
    %define excitatory cell types versus inhibotory types based on trough
    %to peak versus troughdepth/peak height
    
    EI = [alldata(:,5:6)];
    opts = statset('display','final');
    [idx,ctrs] = kmeans(EI,2,'distance','city','Replicates',5,'options',opts);
   
    
   %determine which cells are inhibitory versus which are excitatory
    idx_Duplicate = zeros(size(idx),2);
    idx_Duplicate(:,1)=idx;
    idx_Duplicate(:,2)=idx;
    
%     plot(EI(idx==1,1),EI(idx==1,2),'r.','MarkerSize',12)
%     hold on
%     plot(EI(idx==2,1),EI(idx==2,2),'b.','MarkerSize',12)
%     plot(ctrs(:,1),ctrs(:,2),'kx',...
%         'MarkerSize',12,'LineWidth',2)
%     plot(ctrs(:,1),ctrs(:,2),'ko',...
%         'MarkerSize',12,'LineWidth',2)
%     legend('Cluster 1','Cluster 2','Centroids',...
%         'Location','NW')


    
    for i = 1:size(driftA1,1)
        for j=1:size(driftA1,2)
            driftA1(i,j);
            driftA2(i,j);
            driftB(i,j);
            drift_theta_w(i,j);
            [OSI(i,j) width(i,j) peak(i,j)] = calculate_tuning(driftA1(i,j),driftA2(i,j),driftB(i,j),drift_theta_w(i,j));
            
            
        end
    end
    
   % sprintf('mean firing rate = %f',mean(peak))
 

   
    %variables created to sift through data conditionally 

   responsive = peak(:,1)>2 & peak(:,2)>2;  % firing rate (responsiveness) criteria for whether cells enter subsequent statistical analysis
   responsive_tuned = peak(:,1)>2 & peak(:,2)>2 & OSI(:,1)>0.5 & OSI (:,2)>0.5; % firing rate (responsiveness) criteria for whether cells enter subsequent statistical analysis
   
   
   
   responsive_tuned_E = peak(:,1)>2 & peak(:,2)>2 & OSI(:,1)>0.5 & OSI (:,2)>0.5 & idx_Duplicate(:,1)==2 & idx_Duplicate(:,2)==2; 
   responsive_tuned_I = peak(:,1)>2 & peak(:,2)>2 & OSI(:,1)>0.5 & OSI (:,2)>0.5 & idx_Duplicate(:,1)==1 & idx_Duplicate(:,2)==1;
   
   responsive_tuned_spatial = peak(:,1)>2 & peak(:,2)>2 & driftwpref(:,1) <=0.32 & driftwpref(:,2) <=0.32;
  
   %conditionals to segreagate out excitatory cell pop, I used "idx" as E
   %vs I indentifier, however, it seems to flip between datasets??
   responsive_E = peak(:,1)>2 & peak(:,2)>2 & idx_Duplicate(:,1)==2 & idx_Duplicate(:,2)==2; % responsive units that are excitatory
  
   responsive_E_L2_L3 = peak(:,1)>2 & peak(:,2)>2 & idx_Duplicate(:,1)==2 & idx_Duplicate(:,2)==2 & driftlayer(:,1)<=3 & driftlayer(:,2)<=3;
   %responsive_E_L3 = peak(:,1)>2 & peak(:,2)>2 & idx_Duplicate(:,1)==2 & idx_Duplicate(:,2)==2 & driftlayer(:,1)==3 & driftlayer(:,2)==3;
   responsive_E_L4 = peak(:,1)>2 & peak(:,2)>2 & idx_Duplicate(:,1)==2 & idx_Duplicate(:,2)==2 & driftlayer(:,1)==4 & driftlayer(:,2)==4;
   responsive_E_L5 = peak(:,1)>2 & peak(:,2)>2 & idx_Duplicate(:,1)==2 & idx_Duplicate(:,2)==2 & driftlayer(:,1)==5 & driftlayer(:,2)==5;
   responsive_E_L6 = peak(:,1)>2 & peak(:,2)>2 & idx_Duplicate(:,1)==2 & idx_Duplicate(:,2)==2 & driftlayer(:,1)==6 & driftlayer(:,2)==6;
  
   %conditionals for inhibitory pop
  
   responsive_I = peak(:,1)>2 & peak(:,2)>2 & idx_Duplicate(:,1)==1 & idx_Duplicate(:,2)==1 ; % responsive units that are inhibitory
   responsive_I_L2_L3 = peak(:,1)>2 & peak(:,2)>2 & idx_Duplicate(:,1)==1 & idx_Duplicate(:,2)==1 & driftlayer(:,1)<=3 & driftlayer(:,2)<=3;
   %responsive_I_L3 = peak(:,1)>2 & peak(:,2)>2 & idx_Duplicate(:,1)==1 & idx_Duplicate(:,2)==1 & driftlayer(:,1)==3 & driftlayer(:,2)==3;
   responsive_I_L4 = peak(:,1)>2 & peak(:,2)>2 & idx_Duplicate(:,1)==1 & idx_Duplicate(:,2)==1 & driftlayer(:,1)==4 & driftlayer(:,2)==4;
   responsive_I_L5 = peak(:,1)>2 & peak(:,2)>2 & idx_Duplicate(:,1)==1 & idx_Duplicate(:,2)==1 & driftlayer(:,1)==5 & driftlayer(:,2)==5;
   responsive_I_L6 = peak(:,1)>2 & peak(:,2)>2 & idx_Duplicate(:,1)==1 & idx_Duplicate(:,2)==1 & driftlayer(:,1)==6 & driftlayer(:,2)==6;
   
%    responsive_F1F0 =peak(:,1)>2 & peak(:,2)>2 & driftF1F0(:,1)~=NaN & driftF1F0(:,1)~=NaN;
%    responsive_E_F1F0= peak(:,1)>2 & peak(:,2)>2 & driftF1F0(:,1)~=NaN & driftF1F0(:,1)~=NaN & idx_Duplicate(:,1)==1 & idx_Duplicate(:,2)==1;
%    responsive_I_F1F0=peak(:,1)>2 & peak(:,2)>2 & driftF1F0(:,1)~=NaN & driftF1F0(:,1)~=NaN & idx_Duplicate(:,1)==2 & idx_Duplicate(:,2)==2;
 

      %polar plots for OS at preferred spatial frequency for each data set

   
   
   for moving =1:2 % moving variable describe state and 1:2 describes that there are two states to be considered for this analysis
   meanpeak(dataset,moving) = mean(peak(:,moving));
   meanOSI(dataset,moving) = mean(OSI(responsive,moving));
  
 
   
%    all_OSI=size(OSI);
  
   
   meanOSI_E(dataset,moving) = mean(OSI(responsive_E,moving));
   meanOSI_I(dataset,moving) = mean(OSI(responsive_I,moving));
 
   meanOSI_E_L2_L3(dataset,moving) = mean(OSI(responsive_E_L2_L3,moving));
   %meanOSI_E_L3(dataset,moving) = mean(OSI(responsive_E_L3,moving));
   meanOSI_E_L4(dataset,moving) = mean(OSI(responsive_E_L4,moving));
   meanOSI_E_L5(dataset,moving) = mean(OSI(responsive_E_L5,moving));
   meanOSI_E_L6(dataset,moving) = mean(OSI(responsive_E_L6,moving));
   
   meanOSI_I_L2_L3(dataset,moving) = mean(OSI(responsive_I_L2_L3,moving));
   %meanOSI_I_L3(dataset,moving) = mean(OSI(responsive_I_L3,moving));
   meanOSI_I_L4(dataset,moving) = mean(OSI(responsive_I_L4,moving));
   meanOSI_I_L5(dataset,moving) = mean(OSI(responsive_I_L5,moving));
   meanOSI_I_L6(dataset,moving) = mean(OSI(responsive_I_L6,moving));
   
   
   medianOSI(dataset,moving) = median(OSI(responsive,moving));
   
   meanwidth(dataset,moving) = mean(width(responsive_tuned,moving))
   
   meanspont(dataset,moving) = mean(driftspont(responsive,moving))
   
   meanwpref(dataset,moving) = mean(driftwpref(responsive_tuned_spatial,moving)) 
   
   
   medianwpref(dataset,moving) = median(driftwpref(responsive_tuned_spatial,moving))
   
   meanF1F0(dataset,moving) = mean(driftF1F0(responsive,moving));
   
   %to generate variables to do stats on
   bothwpref{dataset,moving}=driftwpref(responsive,moving);
   bothwpref_E{dataset,moving}=driftwpref(responsive_E,moving);
   bothwpref_I{dataset,moving}=driftwpref(responsive_I,moving);
   
   bothpeak{dataset,moving}=peak(responsive,moving);
   
   bothOSI{dataset,moving}=OSI(responsive,moving);
   bothOSI_E{dataset,moving}=OSI(responsive_E,moving);
   bothOSI_I{dataset,moving} = OSI(responsive_I,moving);
     
   bothOSI_E_L2_L3{dataset,moving}=OSI(responsive_E_L2_L3,moving); 
   bothOSI_E_L4{dataset,moving}=OSI(responsive_E_L4,moving);
   bothOSI_E_L5{dataset,moving}=OSI(responsive_E_L5,moving);
   bothOSI_E_L6{dataset,moving}=OSI(responsive_E_L6,moving);
   
   bothOSI_I_L2_L3{dataset,moving}=OSI(responsive_I_L2_L3,moving);
   bothOSI_I_L4{dataset,moving}=OSI(responsive_I_L4,moving);
   bothOSI_I_L5{dataset,moving}=OSI(responsive_I_L5,moving);
   bothOSI_I_L6{dataset,moving}=OSI(responsive_I_L6,moving);
   
   bothwidth{dataset,moving}=width(responsive_tuned,moving);
   bothwidth_E{dataset,moving}=width(responsive_tuned_E,moving);
   bothwidth_I{dataset,moving}=width(responsive_tuned_I,moving);
   
   bothspont{dataset,moving}=driftspont(responsive,moving);
   
   bothF1F0{dataset,moving}=driftF1F0(responsive,moving);
   bothF1F0_E{dataset,moving}=driftF1F0(responsive_E,moving);
   bothF1F0_I{dataset,moving}=driftF1F0(responsive_I,moving);
   
   
    % to generate error bars for SEM on graphs
    peak_err(dataset,moving) = std(peak(:,moving))/sqrt(length(peak(:,moving)));
    OSI_err(dataset,moving) = std(OSI(responsive,moving))/sqrt(length(OSI(responsive,moving)));
    width_err(dataset,moving) = std(width(responsive_tuned,moving))/sqrt(length(width(responsive_tuned,moving)));
    wpref_err(dataset,moving) = std(driftwpref(responsive_tuned_spatial,moving))/sqrt(length(driftwpref(responsive_tuned_spatial,moving)));
    spont_err(dataset,moving) = std(driftspont(responsive,moving))/sqrt(length(driftspont(responsive,moving)));
    F1F0_err{dataset,moving} = std(driftF1F0(responsive,moving))/sqrt(length(driftF1F0(responsive,moving)));
  end
  
end

%Cris, can you create for loop to generate cumulative distribution plots and variables containing KS test stats for all
%conditions (Adult vs EO1 stationary, seperated out by layers and E/I). want to compare, OSI, DSI, OSI tuning width, SF, F1/F0. Below is my code to generate the plot for one example
%comparison and the code to do stats on the two distributions. Could probably add in the creation of the population polar plots here too

figure

OSI_adult_stationary = bothOSI_E{1};
[f,x,Flo,Fup]= ecdf(OSI_adult_stationary);
plot(x,f,'lineWidth',4,'color','b');
hold on
plot(x,Flo,'b');plot(x,Fup,'b'); %confidence interval for distrib.

hold on

clear f
clear x
clear Flo
clear Fup

OSI_EO1_stationary = bothOSI_I{2};
[f,x,Flo,Fup]= ecdf(OSI_EO1_stationary);
plot(x,f,'lineWidth',4,'color','r');
hold on
plot(x,Flo,'r');plot(x,Fup,'r');
title 'OSI_Excitatory_L2_L3'

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

[pValues_OSI] = ks(OSI_adult_stationary,OSI_EO1_stationary);

% % spatial frequency pref population data comparison
clear f
clear x
clear Flo
clear Fup

figure
wpref_adult_stationary = bothwpref{1};
[f,x,Flo,Fup]= ecdf(wpref_adult_stationary);
plot(x,f,'lineWidth',4,'color','b');
hold on
plot(x,Flo,'b');plot(x,Fup,'b');
hold on

clear f
clear x
clear Flo
clear Fup


wpref_EO1_stationary = bothwpref{2};
[f,x,Flo,Fup]= ecdf(wpref_EO1_stationary);
stairs(x,f,'lineWidth',4,'color','r');
hold on
plot(x,Flo,'r');plot(x,Fup,'r');
hold on
title 'Preferred spatial Frequency'

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

clear f
clear x
clear Flo
clear Fup



%population data on tuning width

figure
tuning_width_adult_stationary = bothwidth{1};
[f,x,Flo,Fup]= ecdf(tuning_width_adult_stationary);
plot(x,f,'lineWidth',4,'color','b');
hold on
plot(x,Flo,'b');plot(x,Fup,'b');
hold on

clear f
clear x
clear Flo
clear Fup

tuning_width_adult_running = bothwidth{3};
[f,x,Flo,Fup]= ecdf(tuning_width_adult_running);
plot(x,f,'lineWidth',4,'color','g');
hold on
plot(x,Flo,'g');plot(x,Fup,'g');
hold on

clear f
clear x
clear Flo
clear Fup

tuning_width_EO1_stationary = bothwidth{2};
[f,x,Flo,Fup]= ecdf(tuning_width_EO1_stationary);
stairs(x,f,'lineWidth',4,'color','r');
hold on
plot(x,Flo,'r');plot(x,Fup,'r');
hold on

clear f
clear x
clear Flo
clear Fup

tuning_width_EO1_running = bothwidth{4};
[f,x,Flo,Fup]= ecdf(tuning_width_EO1_running);
stairs(x,f,'lineWidth',4,'color','m');
hold on
plot(x,Flo,'m');plot(x,Fup,'m');
hold on
title 'OSI tuning width'

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

clear f
clear x
clear Flo
clear Fup

%Bar graphs of the mean
figure
barweb(meanpeak',peak_err');
title('Mean Peak Amplitude')
    %%% print this fig to ps file
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
    %%%%

figure
barweb(meanOSI',OSI_err');
title('Mean OSI')
    %%% print this fig to ps file
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
    %%%%
    
figure
barweb(meanwidth',width_err');
title('Mean Width of Orientation Selective Response')
    %%% print this fig to ps file
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
    %%%%
    
figure
barweb(meanwpref',wpref_err');
title('Mean Preferred Spatial Frequency')
set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

figure
barweb(medianwpref',wpref_err');
title('Median Preferred Spatial Frequency')
set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

figure
barweb(meanspont',spont_err');
title('Mean Spontaneous Activity')
set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

clear h
figure
histbins= 0:2:50;
for rep=1:2
    h(rep,:) = hist(bothspont{rep},histbins)/length(bothspont{rep});
end
bar(histbins,h')
legend('Adult','EO1/2')
title('histogram of spontaneous activity')
ylabel('firing rate') 

%figure
%histbins= -30:2:30;
%for rep=1:2
 %   h(rep,:) = hist(bothpeak{rep})/length(bothpeak{rep});
%end
%bar(histbins,h')
%legend('wt','tg')

clear h
figure
histbins= 0:0.1:1;
for rep=1:2
    h(rep,:) = hist(bothOSI{rep},histbins)/length(bothOSI{rep});
end
bar(histbins,h')
legend('Adult','EO1/2')
title('OSI')

%figure
%histbins= 0:0.2:1;
%for rep=1:2
 %   h(rep,:) = hist(bothwidth{rep},histbins)/length(bothwidth{rep});
%end
%bar(histbins,h')
%legend('wt','tg')


figure
histbins= [0.01 0.02 0.04 0.08 0.16 0.32];
for rep=1:2
    wh(rep,:) = hist(bothwpref{rep},histbins)/length(bothwpref{rep});
end
bar(1:6,wh')
set(gca,'Xtick',1:6);
set(gca,'Xticklabel',histbins)
legend('Adult','EO1/2')
title('spatial frequency preference')

figure
histbins= linspace(0,2,10);
for rep=1:2
    linh(rep,:) = hist(bothF1F0{rep},histbins)/length(bothF1F0{rep});
end
bar(histbins,linh')
legend('Adult','EO1/2')
title('F1F0')

%%% convert ps to pdf and delete old ps
ps2pdf('psfile', psfilename, 'pdffile', [psfilename(1:(end-2)) 'pdf']);
delete(psfilename);


% figure
% barweb(meanF1F0,F1F0_err);
% title('Mean Periodicity')



% Create the area plot using the area function



% figure
% OSI_bins= 0:0.1:1;
% for rep=1:2
%     h(rep,:) = hist(bothOSI{rep},histbins)/length(bothOSI{rep});
% end
% area(OSI_bins,h')
% legend('Adult','EO1/2')
% title('OSI')
% 
% 
% % Add a legend
% legend(groups, 'Location', 'NorthWest');
% 
% % Add title and axis labels
% title('US Population by Age (1860 - 2000)');
% xlabel('Years');
% ylabel('Population in Millions');


ranksum(bothpeak{1},bothpeak{2})
ranksum(bothOSI{1},bothOSI{2})
ranksum(bothOSI_E{1},bothOSI_E{2})
ranksum(bothOSI_I{1},bothOSI_I{2})

% ranksum(bothwpref{1},bothwpref{2})
% ranksum(bothwidth{1},bothwidth{2})
% ranksum(bothspont{1},bothspont{2})
% ranksum(bothF1F0{1},bothF1F0{2})
