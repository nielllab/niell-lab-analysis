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
for dataset = 1:2  %%% adult vs eye opening
    
    
      if dataset ==1
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
            'Good recordings\Adults\11_15_13\rec2\analysis_11_15_13_adult_rec2.mat'}; %% no wn
    elseif dataset ==2
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
            'Good recordings\EO1_EO2\11_27_13_mouseB_EO2\analysis_11_27_13_EO2_mouseB.mat'}; %%% tg
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
        for dontuse =1:1
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
            
            
%                     A1(cellrange,:)=bars_A1;
%                     A2(cellrange,:)=bars_A2;
%                     w(cellrange,:)=bars_w;
%                     theta(cellrange,:)=bars_theta;
%                     bspont(cellrange,:)=barspont;
%             
%                     size(rf_width)
%             
%                     rfw(cellrange,:) = rf_width*30;
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
        
        ev=[]; sp=[];lfp = []; stopcrf=[];mvcrf=[];
        if exist('wn','var')
            for j = 1:length(wn);
             
               if ~isempty(wn(j).N)
                    ev(j) = mean(wn(j).crf(9:12))-mean(wn(j).crf([1 2 19 20]));
                    sp(j) =mean(wn(j).crf([1 2 19 20]));
                    if exist('wn_movement','var')
                        for mv = 1:2
                        lfp(j,mv,:) = interp1(wn_movement(j).freqs,wn_movement(j).mv_lfp(mv,:),1:120);
                        end
                                            stopcrf(j,:)=wn_movement(j).stopCRF;
                   mvcrf(j,:) = wn_movement(j).moveCRF;
                    else
                        lfp=NaN; stopcrf=NaN; mvcrf=NaN;
                        sprintf('no wn movement!!')
                    end

                else
                    
                    ev(j)=NaN;
                    sp(j)=NaN;
                    lfp(j,:,:)=NaN; mvcrf(j,:)=NaN;stopcrf(j,:)=NaN;
                end
            end
            moveLFP(cellrange,:,:)=lfp;
            wn_evoked(cellrange)=ev;
            wn_spont(cellrange)=sp;
            wn_mv(cellrange,:)=mvcrf;
            wn_stop(cellrange,:)=stopcrf;
        else
            display('no wn!!')
            afiles{i}
<<<<<<< HEAD
         wn_mv(cellrange,:)=NaN;
         wn_stop(cellrange,:)=NaN;
=======
            %keyboard
>>>>>>> 160c0cebef96e05de28c5d9f5854ed64fd9595dc
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

EI = [alldata(:,5:6)];
opts = statset('display','final');
[idx,ctrs] = kmeans(EI,2,'distance','city','Replicates',5,'options',opts);
if sum(idx==1)<sum(idx==2)
    inh = (idx==1);
else
    inh = (idx==2);
end



figure
plot(alldata(age==1,5),alldata(age==1,6),'bo'); hold on
plot(alldata(age==2,5),alldata(age==2,6),'go');
legend('EO','adult');

[coeff score latent] = princomp(wvform);
figure
plot(latent);
figure
plot(score(:,1),score(:,2),'o');


inh = alldata(:,6)<1.65 & alldata(:,5)<9;  %%% directly based on wvform; k-means includes another inh group?
midnarrow = alldata(:,6)>1.65 & alldata(:,5)<10;  %%% could analyze these specifically at some point
figure
plot(alldata(find(inh),5),alldata(find(inh),6),'go');
hold on
plot(alldata(find(~inh),5),alldata(find(~inh),6),'bo');

figure
hold on
plot(wvform(find(~inh),:)','b');plot(wvform(find(inh),:)','g');

figure
hold on
plot(score(inh,1),score(inh,2),'go');
plot(score(~inh,1),score(~inh,2),'bo');

figure
hold on
plot(score(midnarrow,1),score(midnarrow,2),'go');
plot(score(~midnarrow,1),score(~midnarrow,2),'bo');

figure
hold on
plot(wvform(age==2 & ~midnarrow',:)','b');
plot(wvform(age==2&midnarrow',:)','g');


inh = (age' ==2 &  alldata(:,6)<1.65 & alldata(:,5)<9) | (age'==1 & wvform(:,16)>0.6);
figure
plot(wvform(age==2,:)','b'); hold on
plot(wvform(age==2 & midnarrow',:)','g');
plot(wvform(age==2 & inh',:)','r');
title('adult waveforms')


figure
plot(wvform(age==1,:)','b'); hold on
plot(wvform(age==1 & midnarrow',:)','g');
plot(wvform(age==1 & inh',:)','r');
title('EO waveforms')

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
<<<<<<< HEAD
keyboard
layerAgePlot(wn_evoked'/30,age,lyr,inh,1,'wn evoked',midnarrow);
layerAgePlot(wn_spont'/30,age,lyr,inh,1,'wn spont',midnarrow);
layerAgeScatter(wn_spont',age,lyr,inh,1,'wn spont',midnarrow);
layerAgeScatter(wn_evoked',age,lyr,inh,1,'wn evoked',midnarrow);
=======
%keyboard
layerAgePlot(wn_evoked/30,age,lyr,inh,1,'wn evoked',midnarrow);
layerAgePlot(wn_spont,age,lyr,inh,1,'wn spont',midnarrow);
layerAgeScatter(wn_spont,age,lyr,inh,1,'wn spont',midnarrow);
layerAgeScatter(wn_evoked,age,lyr,inh,1,'wn evoked',midnarrow);
>>>>>>> 160c0cebef96e05de28c5d9f5854ed64fd9595dc
layerAgeScatter(peak(:,2),age,lyr,inh,1,'moving peak',midnarrow);
figure
hist(wn_evoked); xlabel('wn_evoked')
figure
%%% SF

figure
plot(wn_stop(:,1),wn_mv(:,1),'.'); axis equal; hold on
plot([0 20],[0 20],'g'); xlabel('stationary'); ylabel('moving')
title('spont')

figure
use = ~midnarrow;
plot(wn_stop(use,10)-wn_stop(use,1),wn_mv(use,10)-wn_mv(use,1),'b.'); axis equal;hold on
use = midnarrow;
plot(wn_stop(use,10)-wn_stop(use,1),wn_mv(use,10)-wn_mv(use,1),'r.'); axis equal;hold on
plot([0 20],[0 20],'g');xlabel('stationary'); ylabel('moving')
title('evoked')

layerAgeActivity(wn_stop(:,10)-wn_stop(:,1),wn_mv(:,10)-wn_mv(:,1),age,lyr,inh,1,{'stationary','moving'},'evoked',midnarrow);


LFPmax = max(nanmax(moveLFP,[],3),[],2);
figure
hist(LFPmax)
moveLFP(:,:,59:61)=NaN;



layerAgePlotMv([wn_stop(:,1) wn_mv(:,1)],age,lyr,inh,1,'spont',midnarrow);

%%% peak firing rate in gratings by age and movement)
layerAgeActivity(peak(:,1),peak(:,2),age,lyr,inh,1,{'stationary','moving'},'drift evoked',midnarrow);


figure
hist(peak(:,2),-40:40)
xlabel('drift peak moving');

%%% spont firing rate for gratings by age and movement
layerAgeActivity(driftspont(:,1),driftspont(:,2),age,lyr,inh,1,{'stationary','moving'},'drift spont',midnarrow);


[ m e n] =layerAgePlotMv([wn_stop(:,10)-wn_stop(:,1),wn_mv(:,10)-wn_mv(:,1)],age,lyr,inh,1,'evoked',midnarrow);

[ m e n] =layerAgePlotMv([peak(:,1),peak(:,2)],age,lyr,inh,peak(:,1)>1 & peak(:,2)>1,'drift evoked',midnarrow);

[ m e n] =layerAgePlotMv([driftspont(:,1),driftspont(:,2)],age,lyr,inh,1,'drift spont',midnarrow);

figure
imagesc(squeeze(moveLFP(:,2,:)))

figure
imagesc(squeeze(moveLFP(LFPmax<2*10^4,2,:)))

figure
for ly = 2:6
    plot(1:120,squeeze(nanmedianMW(moveLFP(age'==1 & lyr==ly & LFPmax<2*10^4,1,:),1)),'Linewidth',1+ly/3);
hold on
plot(1:120,squeeze(nanmedianMW(moveLFP(age'==1& lyr==ly& LFPmax<2*10^4,2,:),1)),'g','Linewidth',1+ly/3);
end
title('median EO')

figure
for ly = 2:6
    plot(1:120,squeeze(nanmeanMW(moveLFP(age'==1 & lyr==ly & LFPmax<2*10^4,1,:),1)),'Linewidth',1+ly/3);
hold on
plot(1:120,squeeze(nanmeanMW(moveLFP(age'==1& lyr==ly& LFPmax<2*10^4,2,:),1)),'g','Linewidth',1+ly/3);
end
title('mean EO')


figure
for ly = 2:6
    plot(1:120,squeeze(nanmedianMW(moveLFP(age'==2 & lyr==ly & LFPmax<2*10^4,1,:),1)),'Linewidth',1+ly/3);
hold on
plot(1:120,squeeze(nanmedianMW(moveLFP(age'==2 & lyr==ly & LFPmax<2*10^4,2,:),1)),'g','Linewidth',1+ly/3);
end
title('median adult')

figure
for ly = 2:6
    plot(1:120,squeeze(nanmeanMW(moveLFP(age'==2 & lyr==ly & LFPmax<2*10^4,1,:),1)),'Linewidth',1+ly/3);
hold on
plot(1:120,squeeze(nanmeanMW(moveLFP(age'==2 & lyr==ly & LFPmax<2*10^4,2,:),1)),'g','Linewidth',1+ly/3);
end
title('mean adult')

