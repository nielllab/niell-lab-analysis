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
        
      %  driftF1F0(cellrange,:) = field2array(drift,'F1')./field2array(drift,'F0');
      %  driftF0(cellrange,:) = field2array(drift,'F0');
        %       driftorientfreq_all(cellrange,:)=field2array(drift, 'orientfreq_all');

        driftlayer =  field2array(drift,'layer');
        lyr(cellrange,:) = driftlayer(:,1);
        driftdsi(cellrange,:) = field2array(drift,'dsi');
        %driftOri(cellrange,:) = field2array(drift,'orientfreq_all');
        
%         if exist('bars');
% 
%         bar_spont(cellrange,:)=field2array(bars,'spont');
%        
%         else
%         bar_spont(cellrange,:)= NaN;
%         
%         end
        
        
        
%         if exist('params');
% 
%             
%         all_img_STA(cellrange)= all_img;
%         
%         STA_nx(cellrange)=field2array(params,'nx');
%         STA_ny(cellrange)=field2array(params,'ny');
%         STA_phase(cellrange)=field2array(params,'phase');
%         STA_sigx(cellrange)=field2array(params,'sigx');
%         STA_sigy(cellrange)=field2array(params,'sigy');
%         STA_exp_var(cellrange)=field2array(params,'exp_var');
%         
%         
%         else
%         STA_nx(cellrange)= NaN;
%         STA_ny(cellrange)=  NaN;
%         STA_phase(cellrange)= NaN;
%         STA_sigx(cellrange)= NaN;
%         STA_sigy(cellrange)= NaN; 
%         STA_exp_var(cellrange)=NaN;
% 
%         end
        

        
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
mini= alldata(:,6)>3 & alldata(:,5)>12;
midnarrow = alldata(:,6)>1.75 & alldata(:,5)<10;  %%% could analyze these specifically at some point


figure
plot(alldata(find(inh),5),alldata(find(inh),6),'ro');
hold on
plot(alldata(find(~inh &~mini),5),alldata(find(~inh & ~mini),6),'ko');
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
responsive_stat = peak(:,1)>=1.5;  % firing rate (responsiveness) criteria for whether cells enter subsequent statistical analysis
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
layerAgePlot(peak(:,1),age,lyr,inh,evoked_stat,'drift all evoked Stationary');
peak1_stat=peak(:,1);


%% conduct statistical tests
% normality

[p h]=lillietest(peak(peak(:,1)>=2 & ~inh & age'==2));
g_s=skewness(peak(peak(:,1)>=2 & ~inh & age'==2));
n=size(peak(peak(:,1)>=2 & ~inh & age'==2));

%%to conduct a multiple comparison test of medians bewtween 3 or more
%%groups defined by specific criterion direclty below
age_peak=age(peak(:,1)>=2 & ~inh & lyr>=2);
lyr_peak=lyr(peak(:,1)>=2 & ~inh & lyr>=2);
g=[age_peak',lyr_peak];

g = num2cell(g);
for i=1:size(g,1)
    g{i,1} = [num2str(g{i,1}),num2str(g{i,2})];
end
f={};
f=g(:,1)

 [p,table,stat]=kruskalwallis(peak(peak(:,1)>=2 & ~inh & lyr>=2),f);
 s=multcompare(stat);
 
 %%conduct simple ranksum comparison of medians from two groups directly
 x1=peak(:,1)>=2 & ~inh & lyr==6 & age'==1;
 y1=peak(:,1)>=2 & ~inh & lyr==6 & age'==2;
 
 [p h]=ranksum(peak(x1),peak(y1))
 
 
layerAgePlot_frac_responsive(peak(:,1),age,lyr,inh,responsive_stat,'prct responsive');
%layerAgeScatterMedian(peak1(:,1),age,lyr,inh,evoked_stat,'drift all evoked Stationary');

 %%%%spontaneous rates during drift
% clear driftspont1

driftspont1=driftspont(:,1);
spont_stat_drift = peak(:,1)>=1.5 &  driftspont(:,1)>= 0.045 & driftspont(:,1)<=80;

%layerAgeCDF(driftspont1(:,1),age,lyr,inh,spont_stat_drift,'drift_spont Stationary');
layerAgePlot(driftspont(:,1),age,lyr,inh,spont_stat_drift,'drift_spont Stationary');

clear p h g f g_s age_peak lyr_peak x1 y1
[p h]=lillietest(driftspont1(spont_stat_drift & ~inh & age'==1));
g_s=skewness(driftspont1(spont_stat_drift & ~inh & age'==1))
n=size(driftspont1(spont_stat_drift & ~inh & age'==1))

%%to conduct a multiple comparison test of medians bewtween 3 or more
%%groups defined by specific criterion direclty below
age_peak=age(spont_stat_drift & ~inh & lyr>=2);
lyr_peak=lyr(spont_stat_drift & ~inh & lyr>=2);
g=[age_peak',lyr_peak];

g = num2cell(g);
for i=1:size(g,1)
    g{i,1} = [num2str(g{i,1}),num2str(g{i,2})];
end
f={};
f=g(:,1);

 [p,table,stat]=kruskalwallis(driftspont1(spont_stat_drift & ~inh & lyr>=2),f);
 s=multcompare(stat);
 
 %%conduct simple ranksum comparison of medians from two groups directly
 x1=spont_stat_drift & ~inh & age'==1 & lyr==5;
 y1=spont_stat_drift & ~ inh & age'==2 & lyr==5;
 
 [p h]=ranksum(driftspont1(x1),driftspont1(y1))



tunedOSI_stat =  responsive_stat & OSI (:,1)>=0.5 ; %use=age'==2 & OSI (:,1)>=0.45 & lyr<=3;
tunedOSI_run = responsive_run & OSI (:,2)>=0.5; % firing rate (responsiveness) criteria for whether cells enter subsequent statistical analysis
tunedOSI_either = responsive_either & OSI (:,2)>=0.5|OSI (:,1)>=0.5 ;

%%%prefered orientation
drift_theta_1=size(drift_theta);
drift_theta_1=(drift_theta*180)/pi;
drift_theta_1(drift_theta_1>330)=0;
%drift_theta_1(drift_theta_1>330)=0;

d=drift_theta_1(:,1);
E=d(tunedOSI_stat & age'==1  & lyr<=6 & ~inh)
a=d(tunedOSI_stat & age'==2 & ~inh & lyr<=6)
figure
rose(E,8);
figure
rose(a,8);
figure
hist(a,0:45:330)
figure 
hist(E,0:45:330)


%%preferred oprientation top 30prct
layerAgePlot_pref_Orient(drift_theta_1(:,1),age,lyr,inh,tunedOSI_stat ,{'pref Orient' 'Prct total'},'Prefered Orientation Stationary');


%%%OS

layerAgePlot(OSI(:,1),age,lyr,inh,responsive_stat,'OSI all Stationary');

clear p h g f g_s age_peak lyr_peak x1 y1
[p h]=lillietest(OSI(responsive_stat & ~inh & age'==2));
g_s=skewness(OSI(responsive_stat & ~inh & age'==2))
n=size(OSI(responsive_stat & ~inh & age'==2))

%%to conduct a multiple comparison test of medians bewtween 3 or more
%%groups defined by specific criterion direclty below
age_peak=age(responsive_stat & ~inh & lyr>=2);
lyr_peak=lyr(responsive_stat & ~inh & lyr>=2);
g=[age_peak',lyr_peak];

g = num2cell(g);
for i=1:size(g,1)
    g{i,1} = [num2str(g{i,1}),num2str(g{i,2})];
end
f={};
f=g(:,1);

 [p,table,stat]=kruskalwallis(OSI(responsive_stat & ~inh & lyr>=2),f);
 s=multcompare(stat);
 
 %%conduct simple ranksum comparison of medians from two groups directly
 x1=responsive_stat & ~inh & age'==1 & lyr==6;
 y1=responsive_stat & ~ inh & age'==2 & lyr==6;
 
 [p h]=ranksum(OSI(x1),OSI(y1))
%% calculate unimodal vs multimodal
osi=OSI(responsive_stat & ~inh & lyr>=2 & age'==1);
clear x2
x2 = reshape(osi, 1, prod(size(osi)));
[n, b] = hist(x2, 40);
  
  % This is definitely not probability density function
  x2 = sort(x2);
 % downsampling to speed up computations
 % x2 = interp1 (1:length(x2), x2, 1:1000:length(x2));
  

[dip, p_value, xlow,xup]=HartigansDipSignifTest(x2,500)
%%
% %%%DS
 DSI_stat = driftdsi(:,1)>=0.05  & driftdsi(:,1)<1.05;
 DSI_run =  driftdsi(:,1)>=0.05 & driftdsi(:,2)<1.05;
 tunedDSI_stat = driftdsi(:,1)>=0.5 & driftdsi(:,1)<1.05;
 

 layerAgePlot(driftdsi(:,1),age,lyr,inh,DSI_stat,'DSI Stationary ');

 clear p h g f g_s age_peak lyr_peak x1 y1
 drift=driftdsi(:,1)
[p h]=lillietest(drift(DSI_stat & ~inh & age'==2 & lyr==4));
g_s=skewness(drift(DSI_stat & ~inh & age'==2 & lyr==4))
n=size(drift(DSI_stat & ~inh & age'==2 & lyr==4))

%%to conduct a multiple comparison test of medians bewtween 3 or more
%%groups defined by specific criterion direclty below
age_peak=age(DSI_stat& ~inh & lyr>=3);
lyr_peak=lyr(DSI_stat & ~inh & lyr>=3);
g=[age_peak',lyr_peak];

g = num2cell(g);
for i=1:size(g,1)
    g{i,1} = [num2str(g{i,1}),num2str(g{i,2})];
end
f={};
f=g(:,1);

 [p,table,stat]=kruskalwallis(drift(DSI_stat & ~inh & lyr>=3),f);
 s=multcompare(stat);
 
 %%conduct simple ranksum comparison of medians from two groups directly
 x1=DSI_stat & ~inh & age'==1 & lyr<=3;
 y1=DSI_stat & ~ inh & age'==2 & lyr<=3;
 
 [p h]=ranksum(drift(x1),drift(y1))
 
 
 
 OS=OSI(:,1)>=0.5;
 DS=driftdsi(:,1)>=0.5;
 OS_DS=OSI(:,1)>=0.5 & driftdsi(:,1)>=0.5;

 layerAgePlot_frac_OS_DS(peak(:,1),age,lyr,inh,OS,OS_DS,'OS that is also DS');
 
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
 tunedOSI_stat=OSI(:,1)>=.6 & responsive_stat;
 tunedDSI_stat=driftdsi(:,1)>=0.45 & responsive_stat ;
% use=age'==1  & DSI_tuned; %driftdsi(:,1)>=0.45
 use=age'==1  & lyr==4 & tunedDSI_stat; %driftdsi(:,1)>=0.45
 use=age'==1  & lyr==5 & tunedOSI_stat; %driftdsi(:,1)>=0.45
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
tuned_OS_w_stat =  OSI(:,1)>=0.5 & drift_theta_w(:,1) >=0.05 & drift_theta_w(:,1) <1.3;
layerAgePlot(drift_theta_w(:,1),age,lyr,inh,tuned_OS_w_stat,'OS tuning width Stationary');

OS_w_stat = responsive_stat & drift_theta_w(:,1) >=0.1 & drift_theta_w(:,1) <1.04;
layerAgePlot(drift_theta_w(:,1),age,lyr,inh,OS_w_stat,'OS tuning width Stationary');

drift_theta_w=drift_theta_w(:,1);

clear p h g f g_s age_peak lyr_peak x1 y1
[p h]=lillietest(drift_theta_w(tuned_OS_w_stat & ~inh & age'==2));
g_s=skewness(drift_theta_w(tuned_OS_w_stat& ~inh & age'==2))
n=size(drift_theta_w(tuned_OS_w_stat & ~inh & age'==2))

%%to conduct a multiple comparison test of medians bewtween 3 or more
%%groups defined by specific criterion direclty below
age_peak=age(tuned_OS_w_stat & ~inh & lyr>=2);
lyr_peak=lyr(tuned_OS_w_stat & ~inh & lyr>=2);
g=[age_peak',lyr_peak];

g = num2cell(g);
for i=1:size(g,1)
    g{i,1} = [num2str(g{i,1}),num2str(g{i,2})];
end
f={};
f=g(:,1);

 [p,table,stat]=kruskalwallis(drift_theta_w(tuned_OS_w_stat & ~inh & lyr>=2),f);
 s=multcompare(stat);
 
 %%conduct simple ranksum comparison of medians from two groups directly
 x1=tuned_OS_w_stat& ~inh & age'==1 & lyr==6;
 y1=tuned_OS_w_stat & ~ inh & age'==2 & lyr==6;
 
 [p h]=ranksum(drift_theta_w(x1),drift_theta_w(y1))


%%spatial frequency preference and bandwidth
 
driftwpref(driftwpref==0) = 0.005; %tranform SF pref = 0 to 0.005 in order to plot on log scale
SF_pref_stat = responsive_stat & driftwpref(:,1) >=0.005 & driftwpref(:,1)<=0.40;
%SF_pref_stat_FF= responsive_stat & driftwpref(:,1) ==0.005;
%layerAgeCDF(driftwpref(:,1),age,lyr,inh,SF_pref_stat_either ,'wpref Stationary');
layerAgePlot(driftwpref(:,1),age,lyr,inh,SF_pref_stat,'wpref_stat'); 
layerAgePlot_SF(driftwpref(:,1),age,lyr,inh,SF_pref_stat,'wpref_stat');
driftwpref_stat=driftwpref(:,1);


figure
[f,x]=hist(driftwpref_stat(age'==1 & lyr<=6 & ~inh& SF_pref_stat),0:0.03:0.32);
H1=bar(x,f/sum(f),'b');
title 'EO1 nx L4'
hold on

[f1,x1]=hist(driftwpref_stat(age'==2 & lyr<=6 & ~inh& SF_pref_stat),0:0.03:0.32);
H2=bar(x,f1/sum(f1),'g');
ch=get(H2,'child');
set(ch,'facea',.5)
title 'Eo1_vs adult _ny'
hold on


%%%KS test for difference between distributions
B={driftwpref(SF_pref_stat & age'==1 & lyr<=6),driftwpref(SF_pref_stat & age'==2 & lyr<=6)};
[h p]= kstest2(B{1},B{2})
n1=size(driftwpref(SF_pref_stat & age'==1 & lyr<=6))
n2=size(driftwpref(SF_pref_stat & age'==2 & lyr<=6))

%%dip test for modalness (unimodal or not) EO1
%[x2, n, b] = compute_xpdf(B{:,1});
x_e=B{:,1}
x2 = reshape(x_e, 1, prod(size(x_e)));
  [n, b] = hist(x2, 40); 
 % This is definitely not probability density function
x2 = sort(x2);
[dip, p_value, xlow,xup]=HartigansDipSignifTest(x2,500)

x_A=B{:,2}
x3 = reshape(x_A, 1, prod(size(x_A)));
  [n, b] = hist(x2, 40); 
 % This is definitely not probability density function
x3 = sort(x3);
[dip, p_value, xlow,xup]=HartigansDipSignifTest(x3,500)

clear p h g f g_s age_peak lyr_peak x1 y1
driftwpref_1=driftwpref(:,1);
[p h]=lillietest(driftwpref_1(SF_pref_stat & ~inh & age'==2));
g_s=skewness(driftwpref_1(SF_pref_stat & ~inh & age'==2))
n=size(driftwpref_1(SF_pref_stat & ~inh & age'==2))

%%to conduct a multiple comparison test of medians bewtween 3 or more
%%groups defined by specific criterion direclty below
age_peak=age(SF_pref_stat & ~inh & lyr>=2);
lyr_peak=lyr(SF_pref_stat & ~inh & lyr>=2);
g=[age_peak',lyr_peak];

g = num2cell(g);
for i=1:size(g,1)
    g{i,1} = [num2str(g{i,1}),num2str(g{i,2})];
end
f={};
f=g(:,1);

 [p,table,stat]=kruskalwallis(driftwpref_1(SF_pref_stat & ~inh & lyr>=2),f);
 s=multcompare(stat);
 
 %%conduct simple ranksum comparison of medians from two groups directly
 x1=SF_pref_stat & ~inh & age'==1 & lyr==6;
 y1=SF_pref_stat & ~ inh & age'==2 & lyr==6;
 
 [p h]=ranksum(driftwpref_1(x1),driftwpref_1(y1))


%%%SF_pref_bandwidth

driftwbw_1=driftwbw(:,1);
driftwbw_1(logical(imag(driftwbw_1)))=-1;
%driftwbw_1(driftwbw_1>=7)=0;
driftwbw_1=2*(driftwbw_1);
driftwbw_1(driftwbw_1==-2)=10;
driftwbw_1(driftwbw_1<=0.4)=9;
driftwbw_1(driftwbw_1<=7 & driftwbw_1>=6)=10;

SF_pref_stat =  peak(:,1)>=2 & driftwpref(:,1) >=0.005 & driftwpref(:,1)<=0.40;

clear f x f1 x1
figure
[f,x]=hist(driftwbw_1(age'==1 & SF_pref_stat &lyr>=1 ),0:0.5:10);
H1=bar(x,f/sum(f),'b');
title 'EO1 '
hold on
% y=sum(~isnan(driftwpref_stat(age'==1 & SF_pref_stat_FF)))
% total=y+(sum(~isnan(driftwpref_stat(age'==1 & SF_pref_stat))))
% ff=(y/total)
% bar(ff,x(1,1))
%  hold on
[f1,x1]=hist(driftwbw_1(age'==2 & SF_pref_stat &lyr>=1),0:0.5:10);
H2=bar(x,f1/sum(f1),'g');
ch=get(H2,'child');
set(ch,'facea',.5)
title 'Eo1_vs adult _'

%%test for difference between EO1 and Adult distributions
clear h p B x1 y1
B={driftwbw_1(SF_pref_stat & age'==1 & lyr<=3 & driftwbw_1==9),driftwbw_1(SF_pref_stat & age'==2 & lyr<=3 & driftwbw_1==9)};
[h p]= kstest2(B{1},B{2})
n1= size (B{:,1})
n2= size(B{:,2})

x1=bw_pref & ~inh & age'==1 & lyr<=3;
y1=bw_pref & ~ inh & age'==2 & lyr<=3;
x_size =size(driftwbw_1(bw_pref & ~inh & age'==1 & lyr<=3))
y_size =size(driftwbw_1(bw_pref & ~inh & age'==2 & lyr<=3))
 
 [p h]=ranksum(driftwbw_1(x1),driftwbw_1(y1))
% layerAgeCDF(drift_theta_w(:,2),age,lyr,inh,tuned_OS_w_run ,'OS tuning width Running');

 bw_pref =  peak(:,1)>=2 & driftwpref(:,1) >=0.005 & driftwpref(:,1)<=0.40 & driftwbw_1<6;
SF_pref_stat=peak(:,1)>=2 & driftwpref(:,1) >=0.005 & driftwpref(:,1)<=0.40 ;
layerAgePlot(driftwbw_1(:,1),age,lyr,inh,bw_pref,'SF_pref tuning width Stationary');

% layerAgePlot(drift_theta_w(:,2),age,lyr,inh,tuned_OS_w_run,'OS tuning width Running');
%%% fraction of BW responses that are high pass verus low pass in each
%%% layer
layerAgePlot_SF_BW(driftwbw_1(:,1),age,lyr,inh,SF_pref_stat,'SF_pref tuning width Stationary');
% layerAgeScatterMedian(drift_theta_w(:,2),age,lyr,inh,tuned_OS_w_run,'OS tuning width Running');
j=sum(~isnan(driftwbw_1(peak(:,1)>=2 & driftwbw_1==7)))
clear n b x_e x2 x3 x_A
x_e=B{:,1}
x2 = reshape(x_e, 1, prod(size(x_e)));
  [n, b] = hist(x2, 40); 
 % This is definitely not probability density function
x2 = sort(x2);
 %x2 = interp1 (1:length(x2), x2, 1:1000:length(x2));

[dip, p_value, xlow,xup]=HartigansDipSignifTest(x2,500)

x_A=B{:,2}
x3 = reshape(x_A, 1, prod(size(x_A)));
  [n, b] = hist(x3, 40); 
 % This is definitely not probability density function
x3 = sort(x3);
 %x2 = interp1 (1:length(x2), x2, 1:1000:length(x2));

[dip, p_value, xlow,xup]=HartigansDipSignifTest(x3,500)
%%%Receptive field STA data



%%%Nx data
clear f x f2 x1
good_STA =STA_exp_var>=0.6 &  STA_ny<0.39;


figure
[f,x]=hist(abs(STA_nx(age==1 & lyr'<=6 & ~inh' & good_STA )));
H1=bar(x,f/sum(f),'b');
hold on
[f1,x1]=hist(abs(STA_nx(age==2 & lyr'<=6 & ~inh'& good_STA )));
H2=bar(x,f1/sum(f1),'g');
ch=get(H2,'child');
set(ch,'facea',.5)
title 'Eo1_vs adult _nx'

ranksum(abs(STA_nx(age==1 & lyr'<=6 & ~inh' & good_STA )),abs(STA_nx(age==2 & lyr'<=6 & ~inh'& good_STA )))
v=sum(~isnan(STA_nx(age==1 & lyr'<=6 & ~inh' & good_STA& tunedOSI_stat')))
v1=sum(~isnan(STA_nx(age==2 & lyr'<=6 & ~inh' & good_STA & tunedOSI_stat')))

v2=sum(~isnan(STA_nx(age'==1 & lyr<=6 & ~inh& tunedOSI_stat )))
v3=sum(~isnan(STA_nx(age'==2 & lyr<=6 & ~inh& tunedOSI_stat )))


figure
[f,x]=hist(abs(STA_ny(age==1 & lyr'<=6 & ~inh' & good_STA )));
H1=bar(x,f/sum(f),'b');
hold on
[f1,x1]=hist(abs(STA_ny(age==2 & lyr'<=6 & ~inh'& good_STA )));
H2=bar(x,f1/sum(f1),'g');
ch=get(H2,'child');
set(ch,'facea',.5)
title 'Eo1_vs adult _ny'

clear v v1 v2 v3
v=sum(~isnan(STA_ny(age==1 & lyr'<=6 & ~inh' & good_STA )))
v1=sum(~isnan(STA_ny(age==2 & lyr'<=6 & ~inh' & good_STA )))
                  
%%%layer 4
[f,x]=hist(abs(STA_nx(age==1  & lyr'==4 & good_STA)));
figure
bar(x,f/sum(f));
title 'EO1 nx L4'
 
[f1,x1]=hist(abs(STA_nx(age==2 & lyr'==4 & good_STA)));
figure
bar(x,f1/sum(f1));
title 'adult L4'


%%%plots of ny, nx and ny vs Nx by layer
 jitterValuesX = 2*(rand(size(STA_nx))-0.4)*0.005;   % +/-jitterAmount max
 STA_nx_j=STA_nx+jitterValuesX;

layerAgePlot(abs(STA_nx_j),age,lyr,inh,good_STA','SF_pref tuning width Stationary nx');
layerAgePlot(abs(STA_ny),age,lyr,inh,good_STA','SF_pref tuning width Stationary ny');

STA_nx=STA_nx';
STA_ny=STA_ny';
STA_nx_j=STA_nx_j';
layerAgeNxNy(abs(STA_nx_j),abs(STA_ny),age,lyr,inh,good_STA,{'nx','ny'},'nx vs ny');
layerAgePlot_ratio_jlh(abs(STA_nx),abs(STA_ny),age,lyr,inh,good_STA,{'nx','ny'},'nx vs ny');


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
F1F0_stat = responsive_stat& driftF1F0(:,1) <2.2 & lyr<=6;
F1F0_stat_tuned =tunedOSI_stat & responsive_stat& driftF1F0(:,1) <2.2 & lyr<=6;
% tunedOSI_F1F0_run = tunedOSI_run & driftF1F0(:,2) <=2.1 ;
 %F1F0_run = responsive_run & driftF1F0(:,2) <=2.1;

 
layerAgeCDF(driftF1F0(:,1),age,lyr,inh,F1F0_stat,'F1F0 Stationary');% F1F ratio of all units regardless of selectivity
%layerAgeCDF(driftF1F0(:,2),age,lyr,inh,F1F0_run,'F1F0 Running');
F1F0=driftF1F0(:,1);

clear f x f1 x1

figure
[f,x]=hist(F1F0(age'==1 & F1F0_stat  ),0:0.18:2);
H1=bar(x,f/sum(f),'b');
title 'EO1 '
hold on

[f1,x1]=hist(F1F0(age'==2 & F1F0_stat ),0:0.18:2);
H2=bar(x,f1/sum(f1),'g');
ch=get(H2,'child');
set(ch,'facea',.5);
title 'Eo1_vs adult ';

 
layerAgePlot(driftF1F0(:,1),age,lyr,inh, F1F0_stat,'F1F0 Stationary');
%layerAgePlot(driftF1F0(:,2),age,lyr,inh, F1F0_run,'F1F0 Running');

layerAgePlot_frac_simple(driftF1F0(:,1),age,lyr,inh, F1F0_stat_tuned,'F1F0 Stationary');

%layerAgePlot_frac_simple(driftF1F0(:,2),age,lyr,inh, F1F0_run,'F1F0 Running');

layerAgePlot_frac_simple(driftF1F0(:,1),age,lyr,inh, F1F0_stat,'F1F0 Stationary');
%layerAgePlot_frac_simple(driftF1F0(:,2),age,lyr,inh, F1F0_run,'F1F0 Running');





clear B h p n1 n2 x2 x3 
B={F1F0(F1F0_stat & age'==1 & lyr<=4),F1F0(F1F0_stat & age'==2 & lyr<=4 )};
[h p]= kstest2(B{1},B{2})
n1= size (B{:,1})
n2= size(B{:,2})

x_e=B{:,1}
x2 = reshape(x_e, 1, prod(size(x_e)));
  [n, b] = hist(x2, 40); 
 % This is definitely not probability density function
x2 = sort(x2);
 %x2 = interp1 (1:length(x2), x2, 1:1000:length(x2));

[dip, p_value, xlow,xup]=HartigansDipSignifTest(x2,500)

x_A=B{:,2}
x3 = reshape(x_A, 1, prod(size(x_A)));
  [n, b] = hist(x3, 40); 
 % This is definitely not probability density function
x3 = sort(x3);
 %x2 = interp1 (1:length(x2), x2, 1:1000:length(x2));

[dip, p_value, xlow,xup]=HartigansDipSignifTest(x3,500)

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



