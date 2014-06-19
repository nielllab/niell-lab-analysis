%function compile developmental data
clear all
close all
dbstop if error
% dbclear all % this will exit out of db mode.
%[fname pname] =uiputfile('*.ps','pdf output'); psfilename=fullfile(pname,fname);  %%% get ps filename
psfilename = 'c:/test.ps';   %%% default location
if exist(psfilename,'file')==2;delete(psfilename);end %%% check for previous file

apath = 'D:\Jen_ephys_data\developmental_periods\';
N =0; cells=0;  all_img_STA={};
for dataset = 1:2  %%% adult vs eye opening
    
    if dataset ==1
        
        afiles = {'Good recordings\Adults\1month_2month old\9_25_13\rec1\analysis_9_25_13_P31_rec1.mat',...
            'Good recordings\Adults\1month_2month old\7_18_13_adult\analysis_07_18_13_cluster_rec1.mat',...
            'Good recordings\Adults\5_25_13_Adultmale3mo\Adult B\analysis_5_25_13_B_rec1_adult.mat',...
            'Good recordings\Adults\1month_2month old\4_22_13\analysis_adult.mat',...
            'Good recordings\Adults\7_19_13_adult\analysis_cluster_adult_rec2_7_19_13.mat',...
            'Good recordings\Adults\7_25_13_deep\MouseB\Rec1\analysis_7_25_13_mouseB_rec1.mat',...
            'Good recordings\Adults\7_25_13_deep\MouseB\Rec2\analysis_cluster_data_07_25_13_mouseB_adult_rec2.mat',...
            'Good recordings\Adults\1month_2month old\9_25_13\rec1\analysis_9_25_13_P31_rec1.mat',...
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
%     elseif dataset==2
%         
%         afiles = {'Good recordings\EO7_EO9\8_6_13_EO7\analysis_EO7_rec1_8_6_13.mat',...
%               'Good recordings\EO7_EO9\8_6_13_EO7\rec2\analysis_8_6_13_rec2_EO7.mat',...
%               'Good recordings\EO7_EO9\9_2_13_EO7\rec2\analysis_9_2_13_EO9_rec2.mat',...
%               'Good recordings\EO7_EO9\9_6_13_EO11\analysis_9_6_13_EO11_rec1.mat',...
%               'Good recordings\EO7_EO9\9_6_13_EO11\rec2\analysis_9_6_13_EO11_rec2.mat',...
%               'Good recordings\EO7_EO9\11_27_13_EO7\rec1\analysis_11_27_13_EO7_rec1.mat',...
%               'Good recordings\EO7_EO9\11_27_13_EO7\rec2\analysis_11_27_13_EO7_rec2.mat',...
%               'Good recordings\EO7_EO9\05_10_13_EO7\Record1_upper\analysis_1.mat',...
%               'Good recordings\EO7_EO9\05_10_13_EO7\Record_1_deeper\analysis_EO8_deeper.mat'};
%     elseif dataset ==3
%        
%         afiles = { 'Good recordings\EO3_EO4\5_6_13_EO3\mouseC\analysis_mouseC_EO3_5_6_13.mat',...
%              'Good recordings\EO3_EO4\7_5_13P16_EO4\mouseB\analysis_07_5_13_cluster_7_5_13_mouseB_945uM_analysis.mat',...
%              'Good recordings\EO3_EO4\8_2_13_EO3\analysis_EO3_8_2_13_rec1.mat',...
%              'Good recordings\EO3_EO4\9_10_13_EO3\rec1\analysis_9_10_13_EO3_rec1.mat',...
%              'Good recordings\EO3_EO4\9_10_13_EO3\rec2\analysis_9_10_13_EO3_rec2.mat',...
%              'Good recordings\EO3_EO4\11_25_13_EO4\rec1\analysis_11_25_13_EO4_rec1.mat',...
%              'Good recordings\EO3_EO4\11_25_13_EO4\rec2\analysis_11_25_13_EO4_rec2.mat'}; %% no wn
          
     
    elseif dataset ==2
        afiles = {'Good recordings\EO1_EO2\8_7_13_EO1\rec1_full_clustering\analysis_8_7_13_EO1_rec1_more_strigent.mat',...
            'Good recordings\EO1_EO2\8_7_13_EO1\rec2_full_clustering\analysis_8_7_13_rec2.mat',...
            'Good recordings\EO1_EO2\8_8_13_EO2\rec1_full_clustering\analysis.mat',...
            'Good recordings\EO1_EO2\5_22_13_EO1\analysis_rec1_A_5_22_13_strict_selection.mat',...
            'Good recordings\EO1_EO2\4_29_13_EO1\mouseC\analysis_4_29_13_C.mat',...
            'Good recordings\EO1_EO2\7_17_13_EO1\Analysis_7_17_13_cluster_7_17_13_EO1.mat',...
            'Good recordings\EO1_EO2\7_17_13_EO1\Rec2\Analysis_7_17_13_rec2_EO1.mat',...
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
        
        clear 'rf_width';
        clear params
        clear wn wn_movement
        clear LFP_movement
        clear bars
        
        load([apath afiles{i}]);
        n_units = length(L_ratio);
        cellrange = N+1:N+n_units;
        N=N+n_units;
        
        number(i) = n_units;
        
        alldata( cellrange,1:2) = cells;
        alldata( cellrange,3) = L_ratio;
        
        %%% waveform
        alldata( cellrange,4) = trough_width;
        alldata( cellrange,5) = trough2peak;
        alldata( cellrange,6) = -trough_depth./peak_height;
        alldata( cellrange,7:25)= wv';
        
     
        
        age(cellrange)=3-dataset;
%         for dontuse =1:1
%             %     for c = 1:n_units;
%             %         ch = alldata(c,1);
%             %         cl = alldata(c,2);
%             %         min_t = min(mean_wvform (:,ch : ch+3,cl),[],1);
%             %         [trough_allc trig_chan] = min(min_t);
%             %         trig_chan = ch+trig_chan-1;
%             %         alldata(c,7) = mean_wvform(size(mean_wvform,1),trig_chan,cl)/peak_height(c);
%             %
%             %         t1= squeeze(event_times_all(ch,find(idx_all(ch,:) == cl)));
%             %         dt = diff(t1);
%             %         dt =dt(dt<.02);
%             %         n=hist(dt,.001:0.002:.02);
%             %         [y alldata(c,8)] = max(n);
%             %         n=hist(dt,.0005:0.001:.02);
%             %         alldata(c,9) = max(n(3:8))./mean(n(15:20))
%             %     end;
%             
%             
%                     %   A1(cellrange,:)=bars_A1;
% %                     A2(cellrange,:)=bars_A2;
% %                     w(cellrange,:)=bars_w;
% %                     theta(cellrange,:)=bars_theta;
% %                     bspont(cellrange,:)=barspont;
%                    
%                     
        % end
        
        
        
        if exist('rf_width');
            
        rfw(cellrange,:) = rf_width*30;
        
        else
            
        rfw(cellrange,:) = NaN;
        end
        
%         
        driftA1(cellrange,:)= field2array(drift,'A1');
        driftA2(cellrange,:)=field2array(drift,'A2');
        driftB(cellrange,:)= field2array(drift,'B');
        drift_theta_w(cellrange,:)=field2array(drift,'thetawidth');
        drift_theta(cellrange,:)=field2array(drift,'theta');
        
        driftspont(cellrange,:) = field2array(drift,'spont');
              
        driftwpref(cellrange,:) = field2array(drift,'wpref');
        driftwbw(cellrange,:) = field2array(drift,'bw') ;
        
        driftF1F0(cellrange,:) = field2array(drift,'F1')./field2array(drift,'F0');
        driftF0(cellrange,:) = field2array(drift,'F0');
        %       driftorientfreq_all(cellrange,:)=field2array(drift, 'orientfreq_all');

        driftlayer =  field2array(drift,'layer');
        lyr(cellrange,:) = driftlayer(:,1);
        driftdsi(cellrange,:) = field2array(drift,'dsi');
        %driftOri(cellrange,:) = field2array(drift,'orientfreq_all');
        
        if exist('bars');

        bar_spont(cellrange,:)=field2array(bars,'spont');
       
        else
        bar_spont(cellrange,:)= NaN;
        
        end
        
        
        
        if exist('params');

            
        all_img_STA(cellrange)= all_img;
        
        STA_nx(cellrange)=field2array(params,'nx');
        STA_ny(cellrange)=field2array(params,'ny');
        STA_phase(cellrange)=field2array(params,'phase');
        STA_sigx(cellrange)=field2array(params,'sigx');
        STA_sigy(cellrange)=field2array(params,'sigy');
        STA_exp_var(cellrange)=field2array(params,'exp_var');
        
        
        else
        STA_nx(cellrange)= NaN;
        STA_ny(cellrange)=  NaN;
        STA_phase(cellrange)= NaN;
        STA_sigx(cellrange)= NaN;
        STA_sigy(cellrange)= NaN; 
        STA_exp_var(cellrange)=NaN;

        end
        

        
        %size(wvform)
        wvform(cellrange,:) = wv';
     
        %get firing rate at all measured orients and SF, put into an array:
        %12 rows(orientations) by 7 columns(SpatialFreqs) for each cell
        
        drift_Ori_Sf(cellrange,:) = arrayfun(@(x)(getfield(x,'orientfreq_all')),drift,'UniformOutput',false);
        drift_all(cellrange,:)=drift;
        
       
    end %%% loop over adult vs EO
end
%replace all NaN in F1F0 with 0

% ind = find(isnan(driftF1F0));
% driftF1F0(ind)=0;


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
plot(alldata(age==2,5),alldata(age==2,6),'mo'); hold on


legend('EO','adult');

% [coeff score latent] = princomp(wvform);
% figure
% plot(latent);
% figure
% plot(score(:,1),score(:,2),'o');


inh = alldata(:,6)<1.75 & alldata(:,5)<8.5 ;  %%% directly based on wvform; k-means includes another inh group?
mini= alldata(:,6)>2.5 & alldata(:,5)>12;
midnarrow = alldata(:,6)>1.75 & alldata(:,5)<10;  %%% could analyze these specifically at some point

figure
plot(alldata(find(inh),5),alldata(find(inh),6),'ro');
hold on
plot(alldata(find(~inh),5),alldata(find(~inh),6),'ko');
hold on


figure
plot(alldata(age'==1,5),alldata(age'==1,6),'bo');hold on
plot(alldata(age'==2,5),alldata(age'==2,6),'go')
title 'EO1 vs. adults'



figure
plot(wvform(find(~inh),:)','k');hold on
plot(wvform(find(inh),:)','r'); 


figure
plot(wvform(find(inh),:)','r');hold on
%plot(wvform(age'==2,:)','g');hold on

 y=wvform(find(~inh),:);
 %y=y';
 %plot(1:19,y(2,:));
 
 x=1:19;
 

for i=1:length(y); 
  
    y1= y(i,:);
    
         p  = patchline(x,y1,'edgecolor','k','linewidth',2,'edgealpha',0.1);
end



figure
plot(wvform(age'==1,:)','b');hold on
%plot(wvform(age'==2,:)','g');hold on

 y=wvform(age'==2,:)';
 y=y';
 %plot(1:19,y(2,:));
 
 x=1:19;
 

for i=1:length(y); 
  
    y1= y(i,:);
    
         p  = patchline(x,y1,'edgecolor','g','linewidth',2,'edgealpha',0.1);
end



 


% figure
% hold on
% plot(score(inh,1),score(inh,2),'go');
% plot(score(~inh,1),score(~inh,2),'bo');

% figure
% hold on
% plot(score(midnarrow,1),score(midnarrow,2),'go');
% plot(score(~midnarrow,1),score(~midnarrow,2),'bo');

% figure
% hold on
% plot(wvform(age==2 & ~midnarrow',:)','b');
% plot(wvform(age==2&midnarrow',:)','g');

%Cris's version to hard code a population in the EO group
%as Inhibitory

%inh = (age' ==2 &  alldata(:,6)<1.65 & alldata(:,5)<9) | (age'==1 & wvform(:,16)>0.6);

%inh = (age' ==2 &  alldata(:,6)<1.65 & alldata(:,5)<9);
% 
% figure
% plot(wvform(age==2,:)','g'); hold on
% plot(wvform(age==2 & midnarrow',:)','b');
% plot(wvform(age==2 & inh',:)','r');
% title('adult waveforms')
% 
% 
% figure
% plot(wvform(age==1,:)','g'); hold on
% plot(wvform(age==1 & midnarrow',:)','b');
% plot(wvform(age==1 & inh',:)','r');
% title('EO waveforms')

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

%%firing rate
responsive_run = peak(:,2)>=2; 
responsive_stat = peak(:,1)>=2;  % firing rate (responsiveness) criteria for whether cells enter subsequent statistical analysis
responsive_either= peak(:,2)>=2| peak(:,1)>2;
%responsive_run =  peak(:,2)>=2;
% responsive_both =  peak(:,1)>2 & peak(:,2)>2;



peak1=peak;

peak1(peak1 < 0)=0.01; %%transforms all negative to 0.01 and makes that an equivelent class
%peak1(peak1>0.01 & peak1<1.5 )=1; %%%transforms all 0 values, as "unresponsive" cells to a single class "0.9"
% peak1(peak1==0)=1;
% peak1 = log10(peak1); %%%log base 10 of evoked firing rates

%evoked_either= responsive_either & peak1(:,1)>=0 & peak1(:,1)<=1.6;
evoked_stat= responsive_stat & peak1(:,1)<=80;
evoked_run= responsive_run & peak1(:,2)<=80;

%layerAgeCDF(peak1(:,1),age,lyr,inh,evoked_stat,'drift all evoked Stationary');
layerAgePlot(peak1(:,1),age,lyr,inh,evoked_stat,'drift all evoked Stationary');
peak1_stat=peak(:,1);
age_peak=age(peak(:,1)>=2);
lyr_peak=lyr(peak(:,1)>=2);
g={age_peak';lyr_peak};
p_peak=anovan(peak1_stat(peak(:,1)>=2),g, 'interaction');


layerAgePlot_frac_responsive(peak(:,1),age,lyr,inh,responsive_stat,'prct responsive');
%layerAgeScatterMedian(peak1(:,1),age,lyr,inh,evoked_stat,'drift all evoked Stationary');

 %%%%spontaneous rates during drift
% clear driftspont1

driftspont1=driftspont(:,1);
spont_stat_drift = peak(:,1)>=2 &  driftspont1(:,1)>= 0.01 & driftspont(:,1)<=80;

%layerAgeCDF(driftspont1(:,1),age,lyr,inh,spont_stat_drift,'drift_spont Stationary');
layerAgePlot(driftspont1(:,1),age,lyr,inh,spont_stat_drift,'drift_spont Stationary');


%driftspont1(driftspont1<=.9)=1; %%%transforms zero and close to zero values to 1 for ease of taking log10
%driftspont1=log10(driftspont1);
bar_spont_1=bar_spont(:,1);
resp_stat_spont= responsive_stat & bar_spont_1(:,1)>=0.025; 
%layerAgeScatterMedian(bar_spont_1(:,1),age,lyr,inh,resp_stat_spont,'bar_spont Stationary');
layerAgePlot(bar_spont(:,1),age,lyr,inh,mini,resp_stat_spont,'drift_spont Stationary');
p_peak=anovan(bar_spont_1(peak(:,1)>=2),g, 'interaction');

%ranksum(bar_spont_1(resp_stat_spont & age'==1 & lyr==4),bar_spont_1(resp_stat_spont & age'==2 & lyr==4))
%ranksum(driftspont1(spont_stat_drift & age'==1 & lyr==3),driftspont1(spont_stat_drift & age'==2 & lyr==3))

bar_spont_run=bar_spont(:,2);
resp_run_spont= responsive_run & bar_spont_run(:,1)>=0.045; 
%layerAgeScatterMedian(bar_spont_1(:,1),age,lyr,inh,resp_stat_spont,'bar_spont Stationary');
layerAgePlot(bar_spont_run,age,lyr,inh,mini,resp_run_spont,'bar_spont runnig');
p_peak=anovan(bar_spont_run(peak(:,2)>=2),g, 'interaction');


%%%prefered orientation
drift_theta_1=size(drift_theta);
drift_theta_1=(drift_theta*180)/pi;
drift_theta_1(drift_theta_1>330)=0;
%drift_theta_1(drift_theta_1>330)=0;


% drift_theta_w_1=size(drift_theta_w);
% drift_theta_w_1(drift_theta_w_1<0)=20;
% drift_theta_w_1=(drift_theta_w*180)/pi;

clear top30prct_OSI OSI_stat_top30

top50prct_OSI = prctile(OSI(:,1),50);
OSI_stat_top50 = responsive_either & OSI(:,1)>= top50prct_OSI; 

top50prct_OSI_run = prctile(OSI(:,2),50);
OSI_run_top50 = responsive_run & OSI(:,2)>= top50prct_OSI_run;

tunedOSI_stat =  responsive_stat & OSI (:,1)>=0.5; %use=age'==2 & OSI (:,1)>=0.45 & lyr<=3;
tunedOSI_run = responsive_run & OSI (:,2)>=0.5; % firing rate (responsiveness) criteria for whether cells enter subsequent statistical analysis
tunedOSI_either = responsive_either & OSI (:,2)>=0.5|OSI (:,1)>=0.5 ;

DSI_stat = driftdsi(:,1)>=0 & driftdsi(:,1)<1.05;
DSI_run =  driftdsi(:,2)>=0.05 & driftdsi(:,2)<1.05;
 
 
DSI_tuned= driftdsi(:,1)>=0.5 & driftdsi(:,1)<1.05;
 top50prct_stat = prctile(driftdsi(:,1),50);
 DSI_stat_top50 = responsive_stat & driftdsi(:,1)>= top50prct_stat & driftdsi(:,1)<1.05;
% 
 top50prct_run = prctile(driftdsi(:,2),50);
 DSI_run_top50 = responsive_run & driftdsi(:,2)>= top50prct_run & driftdsi(:,2)<=1;

%%preferred oprientation top 30prct
layerAgePlot_pref_Orient(drift_theta_1(:,1),age,lyr,inh,tunedOSI_stat ,{'pref Orient' 'Prct total'},'Prefered Orientation Stationary');


%%%OS

layerAgePlot(OSI(:,1),age,lyr,inh,DSI_stat,'OSI all Stationary');


%%
% %%%DS
 DSI_stat = driftdsi(:,1)>=0.5 & driftdsi(:,1)<1.05;
 DSI_run =  driftdsi(:,2)>=0.05 & driftdsi(:,2)<1.05;
 
 
DSI_tuned= driftdsi(:,1)>=0.5 & driftdsi(:,1)<1.05;
 top50prct_stat = prctile(driftdsi(:,1),50);
 DSI_stat_top50 = responsive_stat & driftdsi(:,1)>= top50prct_stat & driftdsi(:,1)<1.05; 
 top50prct_run = prctile(driftdsi(:,2),50);
 DSI_run_top50 = responsive_run & driftdsi(:,2)>= top50prct_run & driftdsi(:,2)<=1;
  
% layerAgeCDF(driftdsi(:,1),age,lyr,inh,DSI_stat,'tuned DSI Stationary ');

 layerAgePlot(driftdsi(:,1),age,lyr,inh,mini,DSI_stat,'DSI Stationary ');

 OS=OSI(:,1)>=0.5
 DS=driftdsi(:,1)>=0.5;
 OS_DS=OSI(:,1)>=0.5 & driftdsi(:,1)>=0.5;

 layerAgePlot_frac_OS_DS(driftdsi(:,1),age,lyr,inh,OS,DS,'OS that as DS');
 layerAgePlot(OSI(:,1),age,lyr,inh,mini,DS,'OSI of DS units');

o=OSI(:,1);
d=driftdsi(:,1);
[f,x]=hist(o(age'==1 &  lyr<=5 & DS ),0:0.1:1);
[f1,x1]=hist(o(age'==2 &lyr<=5& DS),0:0.1:1);
bw=[(f/sum(f));(f1/sum(f1))];
figure
bar(x,bw',2)
title 'the OS of DS units'
 
% part 1: create orientations variable
ori = 0:30:330;
ori = ori * pi /180;


%%%polar plot figure
 % use=age'==1 & responsive_stat & lyr==6;
 clear use
% use=age'==1  & DSI_tuned; %driftdsi(:,1)>=0.45
  use=age'==1  & tunedOSI_stat; %driftdsi(:,1)>=0.45

 clear a  ma j a_norm
 
%drift_os=(drift_all(use,1).orientfreq_all);
a = arrayfun(@(x)(getfield(x,'thetatuning')),drift_all(use,1),'UniformOutput',false);

ma=cellfun(@max,a);

figure
for j=1:length(a); 
  
    a_norm = a{j};
    a_norm(a_norm<0)=0;
    a_norm=a_norm./ma(j)';
    if sum(a_norm<0)==0
        polar([ori ori(1)], [a_norm a_norm(1)],'b');hold on
    end
 end

       

%%%OS tuning width
tuned_OS_w_stat =  OSI(:,1)>=0.5 & drift_theta_w(:,1) >=0.1 & drift_theta_w(:,1) <1.2;
%tuned_OS_w_run = tunedOSI_run & drift_theta_w(:,2) >=0.1 & drift_theta_w(:,2) <1.04; % .2 radians =  ~12deg and 1.04radians = ~ 60deg 

%tuned_OS_w_either = tunedOSI_either & drift_theta_w(:,1) >=0.1 & drift_theta_w(:,1) <1.2;

OS_w_stat = OSI_stat_top50 & drift_theta_w(:,1) >=0.1 & drift_theta_w(:,1) <1.04;
%OS_w_run = responsive_run & drift_theta_w(:,2) >=0.1 & drift_theta_w(:,2) <1.04;

%layerAgeCDF(drift_theta_w(:,1),age,lyr,inh,tuned_OS_w_stat ,'OS tuning width Stationary');
% layerAgeCDF(drift_theta_w(:,2),age,lyr,inh,tuned_OS_w_run ,'OS tuning width Running');

layerAgePlot(drift_theta_w(:,1),age,lyr,inh,OS_w_stat,'OS tuning width Stationary');
% layerAgePlot(drift_theta_w(:,2),age,lyr,inh,tuned_OS_w_run,'OS tuning width Running');



%%spatial frequency preference and bandwidth
 
driftwpref(driftwpref==0) = 0.005; %tranform SF pref = 0 to 0.005 in order to plot on log scale


 SF_pref_stat =   driftwpref(:,1) >=0.005 & driftwpref(:,1)<=0.40;

%layerAgeCDF(driftwpref(:,1),age,lyr,inh,SF_pref_stat_either ,'wpref Stationary');


layerAgePlot(driftwpref(:,1),age,lyr,inh,mini,SF_pref_stat,'wpref_stat'); 

driftwpref_stat=driftwpref(:,1);
age_peak=age(peak(:,1)>=2 & SF_pref_stat);
lyr_peak=lyr(peak(:,1)>=2 & SF_pref_stat);
g={age_peak';lyr_peak};
p_peak=anovan(driftwpref_stat(peak(:,1)>=2 &SF_pref_stat),g, 'interaction');


[f,x]=hist(driftwpref_stat(age'==1 & SF_pref_stat & lyr<=4 ),0:0.02:.32)
[f1,x1]=hist(driftwpref_stat(age'==2 & SF_pref_stat & lyr<=4),0:0.02:.32)
bw=[(f/sum(f));(f1/sum(f1))];
figure
bar(x,bw',2)
title 'EO1 vs afult SF pref'

B={driftwpref_stat(SF_pref_stat & age'==1 & lyr<=4),driftwpref_stat(SF_pref_stat & age'==2 & lyr<=4)};

[h p]= kstest2(B{1},B{2});

[f x]= hist(driftwpref_stat(SF_pref_stat& age'==2),0:0.02:0.35);%[f,x]=hist(data(uselist),0:0.02:.35)
figure
bar(x,f/sum(f));

figure
[f x]= hist(driftwpref_stat(SF_pref_stat& age'==1),0:0.02:0.35);%[f,x]=hist(data(uselist),0:0.02:.35)
figure
bar(x,f/sum(f));


layerAgePlot(driftwpref(:,1),age,lyr,inh,mini,SF_pref_stat,'wpref_stat'); 




%%%SF_pref_bandwidth

driftwbw_1=driftwbw(:,1);
driftwbw_1(logical(imag(driftwbw_1)))=-1;
%driftwbw_1(driftwbw_1>=7)=0;
driftwbw_1=2*(driftwbw_1);
driftwbw_1(driftwbw_1==-2)=10;
driftwbw_1(driftwbw_1<=0.4)=9;
driftwbw_1(driftwbw_1<=7 & driftwbw_1>=6)=10;

SF_pref_stat =  peak(:,1)>=0.5 & driftwpref(:,1) >=0.005 & driftwpref(:,1)<=0.32;

[f,x]=hist(driftwbw_1(age'==1 & SF_pref_stat ),0:0.5:10)
[f1,x1]=hist(driftwbw_1(age'==2 & SF_pref_stat ),0:0.5:10)
bw=[(f/sum(f));(f1/sum(f1))];
figure
bar(x,bw',2)
title 'EO1 vs afult bandwidth SF pref'


B={driftwbw_1(SF_pref_stat & age'==1 & lyr==5),driftwbw_1(SF_pref_stat & age'==2 & lyr==5)};

[h p]= kstest2(B{1},B{2});


%layerAgeCDF(driftwbw_1(:,1),age,lyr,inh,SF_pref_stat,'SF_pref tuning width Stationary');
% layerAgeCDF(drift_theta_w(:,2),age,lyr,inh,tuned_OS_w_run ,'OS tuning width Running');
bw_pref =  peak(:,1)>=1 & driftwpref(:,1) >=0.005 & driftwpref(:,1)<=0.32%; & driftwbw_1<7
%bw_pref=SF_pref_stat & driftwbw_1<7;
layerAgePlot(driftwbw_1(:,1),age,lyr,inh,mini,bw_pref,'SF_pref tuning width Stationary');
% layerAgePlot(drift_theta_w(:,2),age,lyr,inh,tuned_OS_w_run,'OS tuning width Running');

%layerAgeScatterMedian(driftwbw_1(:,1),age,lyr,inh,SF_pref_stat,'SF_pref tuning width Stationary');
% layerAgeScatterMedian(drift_theta_w(:,2),age,lyr,inh,tuned_OS_w_run,'OS tuning width Running');

%%%Receptive field STA data

%%%Nx data
clear f x f2 x1
good_STA =STA_exp_var>=0.6;

[f,x]=hist(abs(STA_nx(age==1 & good_STA )),0.08:0.02:0.4);
[f1,x1]=hist(abs(STA_nx(age==2 & good_STA )),0.08:0.02:0.4);
bw=[(f/sum(f));(f1/sum(f1))];
figure
bar(x,bw',2)
title 'EO1 vs afult nx'
v=kstest2(f1,f)

CDF = abs(STA_nx(age==1&good_STA & lyr'==4));
            [f,x,Flo,Fup]= ecdf(CDF);
            stairs(x,f,'lineWidth',4,'color','b');
            hold on
            plot(x,Flo,'color','b');plot(x,Fup,'color','b'); 
         
            hold on
     CDF1 = STA_nx(age==2&good_STA&lyr'==4);
            [f1,x1,Flo1,Fup1]= ecdf(CDF1);
            stairs(x1,f1,'lineWidth',4,'color','g');
            hold on
            plot(x1,Flo1,'color','g');plot(x1,Fup1,'color','g'); 
            axis xy
           
            title 'nx EO1 adult'
                  
%%%all layers
[f,x]=hist(abs(STA_nx(age==1  & lyr'==4 & good_STA)));
figure
bar(x,f/sum(f));
title 'EO1 nx L4'
 
[f1,x1]=hist(abs(STA_nx(age==2 & lyr'==4 & good_STA)));
figure
bar(x,f1/sum(f1));
title 'adult L4'

med_nx_EO1=nanmedian(abs(STA_nx(age==1 & good_STA)))
s_nx_E = semedian(abs(STA_nx(age==1 & good_STA)))

med_nx_adult=nanmedian(abs(STA_nx(age==2  & good_STA)))
s_nx_a = semedian(abs(STA_nx(age==2  & good_STA)))

medians=[med_nx_EO1;med_nx_adult];
SM=[s_nx_E;s_nx_a];

figure
barweb(medians,SM);
title 'median nx'

mean_nx_EO1=nanmean(abs(STA_nx(age==1   & good_STA)))
N=sum(~isnan(STA_nx(age==1  & good_STA)))
SEM_nx_E=nanstd(abs(STA_nx(age==1  & good_STA))/sqrt(N))

mean_nx_adult=nanmean(abs(STA_nx(age==2   & good_STA)))
N=sum(~isnan(STA_nx(age==2  & good_STA)))
SEM_nx_a=nanstd(abs(STA_nx(age==2 & good_STA))/sqrt(N))

means=[mean_nx_EO1;mean_nx_adult];
SEM=[SEM_nx_E;SEM_nx_a];

figure
barweb(means,SM);
title 'means nx'

layerAgePlot(abs(STA_nx),age,lyr,inh,good_STA','SF_pref tuning width Stationary');

%%%plot with fancy multiplot histogram function

% A={(f/sum(f));(f1/sum(f1))}';
% [p h]=kstest2(A{1},A{2});
% ranksum(A{1},A{2})
% %%%all layers
% % A={abs(STA_nx(age==1&good_STA)),abs(STA_nx(age==2&good_STA))};
% % A1={STA_nx(age==1&good_STA),STA_nx(age==2& good_STA)};
% 
% figure
% nhist(A,'legend',{'median=0','median=1'},'median','noerror','separate','smooth');
% title 'all layers'


%%%plots of Ny vs Nx by layer
good_STA=STA_exp_var>=0.6 ;
layerAgeNxNy(abs(STA_nx),abs(STA_ny),age,lyr,inh,good_STA,{'nx','ny'},'nx vs ny');

%%% high SF pref
STA_high =STA_exp_var>=0.6 & driftwpref(:,1)'>0.2;
layerAgeNxNy(abs(STA_nx),abs(STA_ny),age,lyr,inh,STA_high,{'nx','ny'},'nx ratio vs high pref SF(>.10cpd)');

%%% ny nx as a function of spatial frequency pref (regardless of age)
figure
jitterValuesX = 2*(rand(size(STA_nx(lyr'<=6 & STA_high)))-0.5)*.02;   % +/-jitterAmount max
                     %jitterValuesY = 2*(rand(size(data2(uselist)))-0.5)*jitterAmount;   % +/-jitterAmount max
scatter(abs(STA_nx(lyr'<=6 & STA_high))+jitterValuesX, abs(STA_ny(lyr'<=6 & STA_high)), 'ro');hold on
axis square
axis equal
plot([0 0.75],[0 0.75]) 
hold on

STA_low=STA_exp_var>=0.6 & driftwpref(:,1)'<0.02;
jitterValuesX = 2*(rand(size(STA_nx(lyr'<=6 & STA_low)))-0.5)*.02;   % +/-jitterAmount max
                     %jitterValuesY = 2*(rand(size(data2(uselist)))-0.5)*jitterAmount;   % +/-jitterAmount max
scatter(abs(STA_nx(lyr'<=6 & STA_low))+jitterValuesX, abs(STA_ny(lyr'<=6 & STA_low)), 'ko');

x=nanmedian(abs(STA_nx(lyr'<=6 & STA_low)));
x1=nanmedian(abs(STA_nx(lyr'<=6 & STA_high)));

y=nanmedian(abs(STA_ny(lyr'<=6 & STA_low)))
y1=nanmedian(abs(STA_ny(lyr'<=6 & STA_high)))

ranksum(STA_ny(lyr'==4 & STA_low),STA_ny(lyr'==4 & STA_high))

%%%plots of Ny Nx ratios
STA_x_ratio=abs(STA_nx)./abs(STA_ny);
STA_x_ratio(isnan(STA_x_ratio))=0;
good_STA =STA_exp_var>=0.6 & STA_x_ratio>0 ;

%good_STA=good_STA';
layerAgePlot(STA_x_ratio,age,lyr,inh,good_STA,'Nx Ny ratio');


STA_y_ratio=abs(STA_ny)./abs(STA_nx);
good_STA =STA_exp_var>=0.6 %& STA_y_ratio>0 ;
STA_y_ratio(isnan(STA_y_ratio))=0;

good_STA=good_STA';
layerAgePlot(STA_y_ratio,age,lyr,inh,good_STA','Ny Nx ratio');


%%% receptive field size from STAs

%%%all layers

good_STA =STA_exp_var>=0.60;
clear f x f1 x1

% A=(0.5).*STA_sigx;
% B=(0.5).*STA_sigy;
% area=pi.*A.*B;

A= 0.7031.*STA_sigx; %%0.7031 is the calculation for deg per pix 
B= 0.7031.*STA_sigy;

area=pi.*A.*B
figure
hist(area(area<=100))

med_EO1=nanmean(area(age==1   & good_STA))% & area<=100));
Na=sum(~isnan(area(age==1   & good_STA)));
s_EO1= std(area(age==1  & good_STA))/sqrt(Na)% & area<=100 ));


med_adult=nanmedian(area(age==2  & good_STA))% & area<=100 ));
Na1=sum(~isnan(area(age==2   & good_STA)));
s_adult= std(area(age==2   & good_STA))/sqrt(Na1);% & area<=100));

figure
barweb([med_EO1;med_adult],[s_EO1;s_adult]);
title 'median area across all layers'
[p,h]=ranksum(area(age==1  & good_STA & area<=100),area(age==2  & good_STA & area<=100));

% mean_EO1=nanmean(area(age==1 & good_STA ));
% N=sum(~isnan(area(age==1 & good_STA )));
% err_EO1=nanstd(area((age==1 & good_STA )))/sqrt(N);
% 
% mean_adult=nanmean(area(age==2&good_STA ));
% N1=sum(~isnan(area(age==2 & good_STA )));
% err_adult=nanstd(area((age==2 & good_STA )))/sqrt(N1);
% 
% figure
% barweb([mean_EO1;mean_adult],[err_adult;err_EO1]);


x=nanmedian(A(age==1   & good_STA));
sx=semedian(A(age==1 & good_STA));
[p,h]=ranksum(B(age==1  & good_STA ),B(age==2  & good_STA ));

x1=nanmedian(A(age==2  & good_STA));
sx1=semedian(A(age==2 & good_STA));

figure
barweb([x;x1],[sx;sx1]);


% x=nanmean(A(age==1   & good_STA));
% N=sum(~isnan(A(age==1   & good_STA)));
% sx=std(A(age==1 & good_STA))/sqrt(N);
% 
% 
% x1=nanmean(A(age==2  & good_STA));
% N1=sum(~isnan(A(age==2   & good_STA)));
% sx1=std(A(age==2 & good_STA))/sqrt(N1);
% 
% figure
% barweb([x;x1],[sx;sx1])

[p,h]=ranksum((A(age==1 & lyr'==4  & good_STA & area<=100)),(A(age==2 & lyr'==4  & good_STA & area<=100)))




x_B=nanmean(B(age==1   & good_STA));
Nb=sum(~isnan(B(age==1   & good_STA)));
sx_B=std(B(age==1 & good_STA))/sqrt(Nb);

x1_B=nanmean(B(age==2   & good_STA));
Nb2=sum(~isnan(B(age==2   & good_STA)));
sx1_B=std(B(age==2&  good_STA))/sqrt(Nb2);

figure
barweb([x_B;x1_B],[sx_B;sx1_B])

[p,h]=ranksum((B(age==1 & lyr'==4  & good_STA )),(B(age==2 & lyr'==4  & good_STA )))


[f,x]=hist(area(age==1 & good_STA));
figure
bar(x,f/sum(f));
title 'EO1 RF area'

[f,x]=hist(area(age==2 & good_STA));
figure
bar(x,f/sum(f));
title 'adult RF area'


layerAgePlot(area,age,lyr,inh,good_STA','RF area');
layerAgePlot(A,age,lyr,inh,good_STA','RF area');
layerAgePlot(B,age,lyr,inh,good_STA','RF area');



% figure
% hist(B(good_STA))
% 
% figure
% hist(A(good_STA))
% 
% figure
% hist(A(good_STA&age==1));
% 
% figure
% hist(A(good_STA&age==2));




%%%receptive field size(width in radians) from bar stim
% rfw_tuned_stat = tunedOSI_stat & rfw(:,1)<=16& rfw(:,1)>=1;
% %rfw_tuned_run = tunedOSI_run & rfw(:,2)<=16 & rfw(:,2)>=1;
% 
% rfw_all_stat =  responsive_stat & rfw(:,1)<=16 & rfw(:,1)>=1;
% %rfw_all_run =  responsive_run & rfw(:,2)<=16 & rfw(:,2)>=1;
% rfw_all = responsive_either & rfw(:,1)<=16 & rfw(:,1)>=1;
% 
% layerAgeCDF(rfw(:,1),age,lyr,inh,rfw_all_stat,'RFsize Stationary ');
% % layerAgeCDF(rfw(:,2),age,lyr,inh,rfw_all_limit_run,'RFsize Running ');
% 
% layerAgeCDF(rfw(:,1),age,lyr,inh,rfw_tuned_stat,'RFsize Stationary ');
%  
% % layerAgeCDF(rfw(:,1),age,lyr,inh,rfw_tuned_limit_stat,'RFsize Stationary');
% % layerAgeCDF(rfw(:,2),age,lyr,inh,rfw_tuned_limit_run,'RFsize Running ');
% 
% layerAgePlot(rfw(:,1),age,lyr,inh,rfw_tuned_stat,'RFsize Stationary');
% % layerAgePlot(rfw(:,2),age,lyr,inh,rfw_all_limit_run,'RFsize Running ');
% 
% layerAgePlot(rfw(:,1),age,lyr,inh,rfw_all_stat,'RFsize Stationary ');
% % layerAgePlot(rfw(:,1),age,lyr,inh,rfw_tuned_limit_stat,'RFsize Stationary');
% % layerAgePlot(rfw(:,2),age,lyr,inh,rfw_tuned_limit_run,'RFsize Running ');
% 
% layerAgeScatterMedian(rfw(:,1),age,lyr,inh,rfw_all_stat,'RFsize Stationary');
% % layerAgeScatterMedian(rfw(:,2),age,lyr,inh,rfw_all_limit_run,'RFsize Running');
% % 
% layerAgeScatterMedian(rfw(:,1),age,lyr,inh,rfw_tuned_stat,'RFsize Stationary ');
% 
% % layerAgeScatterMedian(rfw(:,1),age,lyr,inh,rfw_tuned_limit_stat,'RFsize Stationary');
% % layerAgeScatterMedian(rfw(:,2),age,lyr,inh,rfw_tuned_limit_run,'RFsize Running ');





%%%simple cells F1F0
tunedOSI_F1F0_stat = tunedOSI_stat & driftF1F0(:,1) <2.2 ;
F1F0_stat = responsive_stat& driftF1F0(:,1) <2.2;

% tunedOSI_F1F0_run = tunedOSI_run & driftF1F0(:,2) <=2.1 ;
 %F1F0_run = responsive_run & driftF1F0(:,2) <=2.1;

 
layerAgeCDF(driftF1F0(:,1),age,lyr,inh,F1F0_stat,'F1F0 Stationary');% F1F ratio of all units regardless of selectivity
%layerAgeCDF(driftF1F0(:,2),age,lyr,inh,F1F0_run,'F1F0 Running');
F1F0=driftF1F0(:,1);

[f,x]=hist(F1F0(age'==1 & ~inh & lyr<5 & F1F0_stat ),0:0.2:2);hold on

[f1,x1]=hist(F1F0(age'==2 & ~inh & lyr<5 &F1F0_stat),0:0.2:2);
bw=[(f/sum(f));(f1/sum(f1))];
%bw=[f;f1]
figure
bar(x,bw',2)
title 'EO1 vs afult F1F0'

[h p]=kstest2(f,f1)

[f, x]=hist(F1F0(F1F0_stat&age'==1 & lyr==4),0:0.1:2);
figure
bar(x,f/sum(f))

[f,x]=hist(F1F0(F1F0_stat&age'==2 & lyr==4),0:0.1:2)
figure
bar(x,f/sum(f))

layerAgeCDF(driftF1F0(:,1),age,lyr,inh,tunedOSI_F1F0_stat,'F1F0 Stationary');
% layerAgeCDF(driftF1F0(:,2),age,lyr,inh,tunedOSI_F1F0_run,'F1F0 Running');
 
layerAgePlot(driftF1F0(:,1),age,lyr,inh,mini, F1F0_stat,'F1F0 Stationary');
%layerAgePlot(driftF1F0(:,2),age,lyr,inh, F1F0_run,'F1F0 Running');

layerAgePlot_frac_simple(driftF1F0(:,1),age,lyr,inh, F1F0_stat,'F1F0 Stationary');
%layerAgePlot_frac_simple(driftF1F0(:,2),age,lyr,inh, F1F0_run,'F1F0 Running');

layerAgePlot_frac_simple(driftF1F0(:,1),age,lyr,inh, tunedOSI_F1F0_stat,'F1F0 Stationary');
%layerAgePlot_frac_simple(driftF1F0(:,2),age,lyr,inh, F1F0_run,'F1F0 Running');



%%% simple cells versus complex cells
layerAgeHist(driftF1F0(:,1),age,lyr,inh, tunedOSI_F1F0_stat,'F1F0 Stationary',midnarrow);
layerAgeHist(driftF1F0(:,2),age,lyr,inh, tunedOSI_F1F0_run,'F1F0 run',midnarrow);

% layerAgeHist(driftwpref(:,1),age,lyr,inh,SFpref_stat,'wpref_Stat',midnarrow);
% layerAgeHist(driftwpref(:,2),age,lyr,inh,SFpref_run,'wpref_run',midnarrow);


figure
hist(lyr(age==1),2:6)
title('EO layer distribution')

figure
hist(lyr(age==2),2:6)
title('adult layer distribution')



%%plot Pref orientation vs. pref SF
driftwpref_1=size(driftwpref);

driftwpref_1=log2(driftwpref);
driftwpref_1=abs(driftwpref_1);



% use = lyr==2|3 & responsive_stat & age'==1 %%& driftwpref_1(:,1)<=3;
% use1= lyr==2|3 & responsive_stat & age'==2 %%& driftwpref_1(:,1)<=3;

use = lyr<=6 & ~inh & tunedOSI_stat & age'==1 %%& driftwpref_1(:,1)<=3;
use1= lyr<=6 & ~inh & tunedOSI_stat & age'==2 %%& driftwpref_1(:,1)<=3;

figure
plot(drift_theta_1(use,1),driftwpref_1(use,1),'bO'); 
axis ij 

figure
plot(drift_theta_1(use1,1),driftwpref_1(use1,1),'gO');hold on
axis ij




%%% peak firing rate in gratings by age and movement)
layerAgeActivity(peak(:,1),peak(:,2),age,lyr,inh,1,{'stationary','moving'},'drift evoked');


figure
hist(peak(:,2),-40:40);hold on
xlabel('drift peak moving');

figure
hist(peak(:,1),-40:40);hold on
xlabel('drift peak stat');

%%% spont firing rate for gratings by age and movement
layerAgeActivity(driftspont(:,1),driftspont(:,2),age,lyr,inh,1,{'stationary','moving'},'drift spont');



