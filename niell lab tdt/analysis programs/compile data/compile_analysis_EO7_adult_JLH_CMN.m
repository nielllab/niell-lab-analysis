%function compile developmental data
clear all
close all
dbstop if error
% dbclear all % this will exit out of db mode.
%[fname pname] =uiputfile('*.ps','pdf output'); psfilename=fullfile(pname,fname);  %%% get ps filename
psfilename = 'c:/test.ps';   %%% default location
if exist(psfilename,'file')==2;delete(psfilename);end %%% check for previous file

apath = 'D:\Jen_ephys_data\developmental_periods\';
N =0; cells=0;
for dataset = 1:2  %%%  EO7 vs adult
    
    
      if dataset ==1
        afiles = {'Good recordings\EO1_EO2\4_29_13_EO1\mouseC\analysis_4_29_13_C.mat',...
            'Good recordings\EO1_EO2\8_7_13_EO1\rec1_full_clustering\analysis_8_7_13_EO1_rec1_more_strigent.mat',...
            'Good recordings\EO1_EO2\8_7_13_EO1\rec2_full_clustering\analysis_8_7_13_rec2.mat',...
            'Good recordings\EO1_EO2\5_22_13_EO1\analysis_rec1_A_5_22_13_strict_selection.mat',...
            'Good recordings\EO1_EO2\7_17_13_EO1\Analysis_7_17_13_cluster_7_17_13_EO1.mat',...
            'Good recordings\EO1_EO2\7_17_13_EO1\Rec2\Analysis_7_17_13_rec2_EO1.mat',...
            'Good recordings\EO1_EO2\8_8_13_EO2\rec1_full_clustering\analysis.mat',...
            'Good recordings\EO1_EO2\9_9_13_EO1\rec1\analysis_9_9_13_EO1_rec1.mat',...
            'Good recordings\EO1_EO2\9_30_13_EO1\rec1\analysis_9_30_13_EO1_rec1.mat',...
            'Good recordings\EO1_EO2\9_30_13_EO1\rec2\analysis_9_30_13_EO1_rec2.mat',...
            'Good recordings\EO1_EO2\9_9_13_EO1\rec2\analysis_9_9_13_EO1_rec2.mat',...
            'Good recordings\EO1_EO2\10_1_13_EO2\rec1\analysis_10_1_13_EO2_rec1.mat',...
            'Good recordings\EO1_EO2\11_20_13_EO0\rec1\analysis_11_20_13_rec1.mat',...
            'Good recordings\EO1_EO2\11_20_13_EO0\rec2\analysis_11_20_13_EO0_rec2.mat',...
            'Good recordings\EO1_EO2\11_21_13_EO1_mouseA\analysis_11_21_13_mouseA_EO1.mat',...
            'Good recordings\EO1_EO2\11_23_13_EO2\analysis_11_23_13_EO2.mat',...
            'Good recordings\EO1_EO2\11_26_13_mouseA_EO1\rec1\analysis_11_26_13_EO1_mouseA_rec1.mat',...
            'Good recordings\EO1_EO2\11_26_13_mouseA_EO1\rec2\analysis_11_26_13_EO1_rec2.mat',...
            'Good recordings\EO1_EO2\11_26_13_mouseB_EO1\analysis_11_26_13_mouseB_EO1.mat',...
            'Good recordings\EO1_EO2\11_27_13_mouseB_EO2\analysis_11_27_13_EO2_mouseB.mat'}; %% no wn
    elseif dataset ==2
        afiles = {'Good recordings\Adults\1month_2month old\4_22_13\analysis_adult.mat',...
            'Good recordings\Adults\5_25_13_Adultmale3mo\Adult B\analysis_5_25_13_B_rec1_adult.mat',...
            'Good recordings\Adults\1month_2month old\7_18_13_adult\analysis_07_18_13_cluster_rec1.mat',...
            'Good recordings\Adults\7_19_13_adult\analysis_cluster_adult_rec2_7_19_13.mat',...
            'Good recordings\Adults\7_25_13_deep\MouseB\Rec1\analysis_7_25_13_mouseB_rec1.mat',...
            'Good recordings\Adults\7_25_13_deep\MouseB\Rec2\analysis_cluster_data_07_25_13_mouseB_adult_rec2.mat',...
            'Good recordings\Adults\1month_2month old\9_25_13\rec1\analysis_9_25_13_P31_rec1.mat',...
            'Good recordings\Adults\9_27_13\rec1\analysis_9_27_13_adult_rec1.mat',...
            'Good recordings\Adults\9_27_13\rec2\analysis_9_27_13_adult_rec2.mat',...
            'Good recordings\Adults\9_12_13\rec1\analysis_9_12_13_adult_rec1.mat',...
            'Good recordings\Adults\9_12_13\rec2\analysis_9_12_13_adult_rec2.mat',...
            'Good recordings\Adults\11_11_13\rec1\analysis_11_11_13_adult_rec1.mat',... %%% no wn
            'Good recordings\Adults\11_11_13\rec2\analysis_11_11_13_adult_rec2.mat',...
            'Good recordings\Adults\11_13_13\rec1\analysis_adult_11_13_13_rec1.mat',...
            'Good recordings\Adults\11_13_13\rec2\analysis_11_13_13_rec2.mat',...
            'Good recordings\Adults\11_14_13\rec1\analysis_11_14_13_adult_rec1.mat',...
            'Good recordings\Adults\11_14_13\rec2\analysis_11_14_13_adult_rec2.mat',...
            'Good recordings\Adults\11_15_13\rec1\analysis_11_15_13_adult_rec1.mat',...
            'Good recordings\Adults\11_15_13\rec2\analysis_11_15_13_adult_rec2.mat'}; %%% tg
        
     
          
    end
    
    
    for i = 1:length(afiles)
        

        clear wn wn_movement
        
        load([apath afiles{i}]);
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
        
        age(cellrange)=3-dataset;
        for dontuse =1:0
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
        end
        
        driftA1(cellrange,:)= field2array(drift,'A1');
        driftA2(cellrange,:)=field2array(drift,'A2');
        driftB(cellrange,:)= field2array(drift,'B');
        drift_theta_w(cellrange,:)=field2array(drift,'thetawidth');
        driftspont(cellrange,:) = field2array(drift,'spont');
              
        driftwpref(cellrange,:) = field2array(drift,'wpref');
        driftwbw(cellrange,:) = field2array(drift,'bw') ;
        
        driftF1F0(cellrange,:) = field2array(drift,'F1')./field2array(drift,'F0');
        driftF0(cellrange,:) = field2array(drift,'F0');
        %       driftorientfreq_all(cellrange,:)=field2array(drift, 'orientfreq_all');
        
        driftlayer =  field2array(drift,'layer');
        lyr(cellrange,:) = driftlayer(:,1);
        driftdsi(cellrange,:) = field2array(drift,'dsi');
        
        %size(wvform)
        wvform(cellrange,:) = wv';
        
        ev=[]; sp=[];lfp = [];
        if exist('wn','var')
            for j = 1:length(wn);
             
               if ~isempty(wn(j).N)
                    ev(j) = mean(wn(j).crf(9:12))-mean(wn(j).crf([1 2 19 20]));
                    sp(j) =mean(wn(j).crf([1 2 19 20]));
                    if exist('wn_movement','var')
                        for mv = 1:2
                        lfp(j,mv,:) = interp1(wn_movement(j).freqs,wn_movement(j).mv_lfp(mv,:),1:120);
                        end
                    else
                        lfp=NaN;
                        sprintf('no wn movement!!')
                    end
%                     stopcrf(i,:)=wn_movement(i).stopCRF;
%                     mvcrf(i,:) = wn_movement(i).moveCRF;
                else
                    
                    ev(j)=NaN;
                    sp(j)=NaN;
                    lfp(j,:,:)=NaN;
                end
            end
            moveLFP(cellrange,:,:)=lfp;
            wn_evoked(cellrange)=ev;
            wn_spont(cellrange)=sp;
        else
            display('no wn!!')
            afiles{i}
            %keyboard
           moveLFP(cellrange,:,:)=NaN;
           wn_evoked(cellrange)=NaN;
            wn_spont(cellrange)=NaN;
        end
        size(lfp);
      
    end   %%%  loop over afiles
end
%%% loop over adult vs EO


%replace all NaN in F1F0 with 0

ind = find(isnan(driftF1F0));
driftF1F0(ind)=0;


%define excitatory cell types versus inhibotory types based on trough
%to peak versus troughdepth/peak height

EI = [alldata(:,4:6)];
opts = statset('display','final');
[idx,ctrs] = kmeans(EI,3,'distance','city','Replicates',5,'options',opts);
% if sum(idx==1)<sum(idx==2)
%     inh = (idx==1);
% else
%     inh = (idx==2);
% end

figure
plot(alldata(:,5),alldata(:,6),'ko') %all waveforms included in the analysis are clustered here
axis([0 15 0 5])

%plot scatter of all cells seperated by three waveform characteristics:p to
%trough, p:t height ratio and width of trough

figure 
scatter3(alldata(idx==1,5),alldata(idx==1,6),alldata(idx==1,4),'bo'); hold on
scatter3(alldata(idx==2,5),alldata(idx==2,6),alldata(idx==2,4),'go'); hold on
scatter3(alldata(idx==3,5),alldata(idx==3,6),alldata(idx==3,4),'ro');
axis([0 15 0 5 0 6])




inh = (idx==3);
midnarrow = (idx==1);
broad = (idx==2);

% figure
% scatter3(alldata(age==1,5),alldata(age==1& idx==1,6),alldata(age==1& idx==1,4),'bo'); hold on
% scatter3(alldata(age==2,5),alldata(age==2,6),alldata(age==2,4),'rx');
% legend('EO','adult');

[coeff score latent] = princomp(wvform);
figure
plot(latent);
figure
plot(score(:,3),score(:,4),'o');


% inh = alldata(:,6)<1.65 & alldata(:,5)<9;  %%% directly based on wvform; k-means includes another inh group?
% midnarrow = alldata(:,6)>1.65 & alldata(:,5)<10;  %%% could analyze these specifically at some point
figure
plot(alldata(find(inh),5),alldata(find(inh),6),'ro');
hold on
plot(alldata(find(broad),5),alldata(find(broad),6),'go'); hold on
%plot(alldata(find(midnarrow),5),alldata(find(midnarrow),6),'bo');
axis ([0 15 0 5])


% figure
% plot(alldata(find(midnarrow),5),alldata(find(midnarrow),6),'bo');
% ylabel('Peak:trough ratio')
% xlabel('trough to peak distance')

figure
hold on
plot(wvform(find(inh),:)','r');plot(wvform(find(broad),:)','g');
hold on
plot(wvform(find(midnarrow),:)','b')
ylabel('amplitude')
xlabel('msec')
% 
% figure
% hold on
% plot(wvform(find(broad),:)','g');plot(wvform(find(inh),:)','r'); hold on
% plot(wvform(find(midnarrow),:)','b')
% ylabel('amplitude')
% xlabel('msec')
% 
% figure
% hold on
% plot(score(inh,1),score(inh,2),'go');
% plot(score(~inh,1),score(~inh,2),'bo');
% 
% figure
% hold on
% plot(score(midnarrow,1),score(midnarrow,2),'go');
% plot(score(~midnarrow,1),score(~midnarrow,2),'bo');

figure
hold on
plot(wvform(age==1 & inh',:)','g');
hold on
plot(wvform(age==2 & inh',:)','b');%age==2 = EO7 and age==1 = adult plots all inhibotry waveforms as a function of age


figure
hold on
plot(wvform(age==1 & midnarrow',:)','g');
hold on
plot(wvform(age==2 & midnarrow',:)','b');

figure
hold on
plot(wvform(age==1 & broad',:)','g');
hold on
plot(wvform(age==2 & broad',:)','b');


% inh = (age' ==2 &  alldata(:,6)<1.65 & alldata(:,5)<9) | (age'==1 & wvform(:,16)>0.6);
% figure
% plot(wvform(age==2,:)','b'); hold on
% %plot(wvform(age==2 & midnarrow',:)','g');
% plot(wvform(age==2 & inh',:)','r');
% title('EO3 waveforms')


figure
plot(wvform(age==2,:)','b'); hold on
%plot(wvform(age==2 & midnarrow',:)','g');
plot(wvform(age==2 & inh',:)','r');
title('EO3 waveforms')

figure
plot(wvform(age==1,:)','b'); hold on
%plot(wvform(age==1 & midnarrow',:)','g');
plot(wvform(age==1 & inh',:)','r');
title('EO7 waveforms')

%%% spot analysis, e.g. size, on/off, latency

%%% STA analysis - svd, fit to gabors/gaussian = which fits better?


for i = 1:size(driftA1,1)
    for j=1:size(driftA1,2)
        driftA1(i,j);
        driftA2(i,j);
        driftB(i,j);
        drift_theta_w(i,j);
        [OSI(i,j) width(i,j) peak(i,j)] = calculate_tuning(driftA1(i,j),driftA2(i,j),driftB(i,j),drift_theta_w(i,j));
        
        
    end
end

%variables created to sift through data conditionally

responsive = peak(:,1)>2 & peak(:,2)>2;  % firing rate (responsiveness) criteria for whether cells enter subsequent statistical analysis
tunedOSI =OSI(:,1)>0.5 & OSI (:,2)>0.5; % firing rate (responsiveness) criteria for whether cells enter subsequent statistical analysis

tunedSF = driftwpref(:,1) <=0.32 & driftwpref(:,2) <=0.32;


%Cris, can you create for loop to generate cumulative distribution plots and variables containing KS test stats for all
%conditions (Adult vs EO1 stationary, seperated out by layers and E/I). want to compare, OSI, DSI, OSI tuning width, SF, F1/F0. Below is my code to generate the plot for one example
%comparison and the code to do stats on the two distributions. Could probably add in the creation of the population polar plots here too

figure
hist(lyr(age==1),2:6)
title('EO layer distribution')

figure
hist(lyr(age==2),2:6)
title('adult layer distribution')



layerAgePlot(OSI(:,2),age,lyr,inh,responsive,'OSI',midnarrow);
layerAgePlot(driftwpref(:,2),age,lyr,inh,responsive,'wpref',midnarrow);
layerAgeScatter(driftwpref(:,2),age,lyr,inh,responsive,'wpref',midnarrow);
%keyboard
layerAgePlot(wn_evoked/30,age,lyr,inh,1,'wn evoked',midnarrow);
layerAgePlot(wn_spont,age,lyr,inh,1,'wn spont',midnarrow);
layerAgeScatter(wn_spont,age,lyr,inh,1,'wn spont',midnarrow);
layerAgeScatter(wn_evoked,age,lyr,inh,1,'wn evoked',midnarrow);
layerAgeScatter(peak(:,2),age,lyr,inh,1,'moving peak',midnarrow);
figure
hist(wn_evoked); xlabel('wn_evoked')
figure
%%% SF

%%% fix normalization
figure
for ly = 2:5
    plot(1:120,squeeze(nanmean(moveLFP(age'==1 & lyr==ly,1,:),1)));
hold on
plot(1:120,squeeze(nanmean(moveLFP(age'==1& lyr==ly,2,:),1)),'g');
end

figure
for ly = 2:5
    plot(1:120,squeeze(nanmean(moveLFP(age'==2 & lyr==ly,1,:),1)));
hold on
plot(1:120,squeeze(nanmean(moveLFP(age'==2 & lyr==ly,2,:),1)),'g');
end

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

OSI_EO1_stationary = bothOSI_E{2}; %%% bug here - was bothOSI_I{2}
[f,x,Flo,Fup]= ecdf(OSI_EO1_stationary);
plot(x,f,'lineWidth',4,'color','r');
hold on
plot(x,Flo,'r');plot(x,Fup,'r');
title 'OSI_Excitatory_L2_L3'

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

[pValues_OSI] = ks(OSI_adult_stationary,OSI_EO1_stationary);

% % spatial frequency pref population data comparison

figure
wpref_adult_stationary = bothwpref{1};
[f,x,Flo,Fup]= ecdf(wpref_adult_stationary);
plot(x,f,'lineWidth',4,'color','b');
hold on
plot(x,Flo,'b');plot(x,Fup,'b');
hold on

wpref_EO1_stationary = bothwpref{2};
[f,x,Flo,Fup]= ecdf(wpref_EO1_stationary);
stairs(x,f,'lineWidth',4,'color','r');
hold on
plot(x,Flo,'r');plot(x,Fup,'r');
hold on
title 'Preferred spatial Frequency'

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');



%population data on tuning width

figure
tuning_width_adult_stationary = bothwidth{1};
[f,x,Flo,Fup]= ecdf(tuning_width_adult_stationary);
plot(x,f,'lineWidth',4,'color','b');
hold on
plot(x,Flo,'b');plot(x,Fup,'b');
hold on

tuning_width_adult_running = bothwidth{3};
[f,x,Flo,Fup]= ecdf(tuning_width_adult_running);
plot(x,f,'lineWidth',4,'color','g');
hold on
plot(x,Flo,'g');plot(x,Fup,'g');
hold on

tuning_width_EO1_stationary = bothwidth{2};
[f,x,Flo,Fup]= ecdf(tuning_width_EO1_stationary);
stairs(x,f,'lineWidth',4,'color','r');
hold on
plot(x,Flo,'r');plot(x,Fup,'r');
hold on

tuning_width_EO1_running = bothwidth{4};
[f,x,Flo,Fup]= ecdf(tuning_width_EO1_running);
stairs(x,f,'lineWidth',4,'color','m');
hold on
plot(x,Flo,'m');plot(x,Fup,'m');
hold on
title 'OSI tuning width'

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');


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
