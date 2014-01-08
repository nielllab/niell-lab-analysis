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
        afiles = {'Good recordings\Adults\5_25_13_Adultmale3mo\Adult B\analysis_5_25_13_B_rec1_adult.mat',...
            'Good recordings\Adults\1month_2month old\7_18_13_adult\analysis_07_18_13_cluster_rec1.mat',...
            'Good recordings\Adults\1month_2month old\4_22_13\analysis_adult.mat',...
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
        afiles = {'Good recordings\EO1_EO2\8_7_13_EO1\rec1_full_clustering\analysis_8_7_13_EO1_rec1_more_strigent.mat',...
            'Good recordings\EO1_EO2\8_7_13_EO1\rec2_full_clustering\analysis_8_7_13_rec2.mat',...
            'Good recordings\EO1_EO2\5_22_13_EO1\analysis_rec1_A_5_22_13_strict_selection.mat',...
            'Good recordings\EO1_EO2\8_8_13_EO2\rec1_full_clustering\analysis.mat',...
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
        

        clear wn wn_movement
        clear LFP_movement
        
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
            
            
                    %   A1(cellrange,:)=bars_A1;
%                     A2(cellrange,:)=bars_A2;
%                     w(cellrange,:)=bars_w;
%                     theta(cellrange,:)=bars_theta;
%                     bspont(cellrange,:)=barspont;
                   
                    
        end
        
        %size(rf_width)
        
        if exist('rf_width');
            
        rfw(cellrange,:) = rf_width*30;
        
        else
            
        rfw(cellrange,:) = NaN;
        end
        clear 'rf_width';
        
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

        %size(wvform)
        wvform(cellrange,:) = wv';
     
        %get firing rate at all measured orients and SF, put into an array:
        %12 rows(orientations) by 7 columns(SpatialFreqs) for each cell
        
        drift_Ori_Sf(cellrange,:) = arrayfun(@(x)(getfield(x,'orientfreq_all')),drift,'UniformOutput',false);
        drift_all(cellrange,:)=drift;
        
        ev=[]; sp=[];lfp = []; stopcrf=[];mvcrf=[];lfp_bar = [];
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
% <<<<<<< HEAD
         wn_mv(cellrange,:)=NaN;
         wn_stop(cellrange,:)=NaN;
% =======
%             %keyboard
% >>>>>>> 160c0cebef96e05de28c5d9f5854ed64fd9595dc
           moveLFP(cellrange,:,:)=NaN;
           wn_evoked(cellrange)=NaN;
           wn_spont(cellrange)=NaN;
         
         
        end
          
   if exist ('LFP_movement', 'var')
    
      for k = 1:length(LFP_movement)   
           if ~isempty(LFP_movement(k).mv_lfp) 
           for mv = 1:2
                lfp_bar(k,mv,:) = interp1(LFP_movement(k).freqs,LFP_movement(k).mv_lfp(mv,:),1:120);
           end 
           else 
                lfp_bar(k,:,:) = NaN;
                sprintf('no LFP_movement');
           end     
      end
  
   moveLFP_bar(cellrange,:,:)=lfp_bar;  
   
  else
       
   moveLFP_bar(cellrange,:,:)=NaN;
  
   end
  
    end %%loop over cells
         
end %%% loop over adult vs EO

           
    


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
plot(alldata(age==2,5),alldata(age==2,6),'go');
legend('EO','adult');

% [coeff score latent] = princomp(wvform);
% figure
% plot(latent);
% figure
% plot(score(:,1),score(:,2),'o');


inh = alldata(:,6)<1.75 & alldata(:,5)<8.5;  %%% directly based on wvform; k-means includes another inh group?
midnarrow = alldata(:,6)>1.75 & alldata(:,5)<10;  %%% could analyze these specifically at some point

figure
plot(alldata(find(inh),5),alldata(find(inh),6),'ro');
hold on
plot(alldata(find(~inh),5),alldata(find(~inh),6),'ko');
hold on

figure
plot(alldata(age'==2,5),alldata(age'==2,6),'go');hold on
plot(alldata(age'==1,5),alldata(age'==1,6),'bo')


figure
hold on
plot(wvform(find(~inh),:)','k');plot(wvform(find(inh),:)','r');

figure
plot(wvform(age'==2,:)','g');hold on
plot(wvform(age'==1,:)','b');
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

inh = (age' ==2 &  alldata(:,6)<1.65 & alldata(:,5)<9);
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
responsive_either = peak(:,1)>=2 | peak(:,2)>=2;
responsive_stat = peak(:,1)>2;  % firing rate (responsiveness) criteria for whether cells enter subsequent statistical analysis
%responsive_run =  peak(:,2)>=2;
% responsive_both =  peak(:,1)>2 & peak(:,2)>2;


peak1 = size(peak);
peak1=peak;

peak1(peak1 < 0)=0.01; %%transforms all negative to 0.01 and makes that an equivelent class
peak1(peak1>0.01 & peak1<1.5 )=1; %%%transforms all 0 values, as "unresponsive" cells to a single class "0.9"
peak1(peak1==0)=1;
peak1 = log10(peak1); %%%log base 10 of evoked firing rates

%evoked_either= responsive_either & peak1(:,1)>=0 & peak1(:,1)<=1.6;
evoked_stat= responsive_stat & peak1 (:,1)>=0 & peak1(:,1)<=40;

layerAgeCDF(peak1(:,1),age,lyr,inh,evoked_stat,'drift all evoked Stationary');
layerAgePlot(peak1(:,1),age,lyr,inh,evoked_stat,'drift all evoked Stationary');
layerAgeScatterMedian(peak1(:,1),age,lyr,inh,evoked_stat,'drift all evoked Stationary');

 %%%%spontaneous rates during drift
% clear driftspont1
driftspont1 = size(driftspont);
driftspont1=driftspont;
driftspont1(driftspont1<=.9)=1; %%%transforms zero and close to zero values to 1 for ease of taking log10
driftspont1=log10(driftspont1);

spont_stat = peak(:,1)>=2 &  driftspont1(:,1)>=0 &  driftspont1(:,1)<=20 ;

layerAgeCDF(driftspont1(:,1),age,lyr,inh,spont_stat,'drift_spont Stationary');
layerAgePlot(driftspont1(:,1),age,lyr,inh,spont_stat,'drift_spont Stationary');
layerAgeScatterMedian(driftspont1(:,1),age,lyr,inh,spont_stat,'drift_spont Stationary');



%%%prefered orientation
drift_theta_1=size(drift_theta);
drift_theta_1=(drift_theta*180)/pi;
drift_theta_w_1=size(drift_theta_w);
drift_theta_w_1(drift_theta_w_1<0)=20;

drift_theta_w_1=(drift_theta_w*180)/pi;

layerAgeCDF(drift_theta_1(:,1),age,lyr,inh,responsive_stat,'drift all evoked Stationary');
layerAgeScatterMedian(drift_theta_1(:,1),age,lyr,inh,responsive_stat,'drift all evoked Stationary');
layerAgePlot(drift_theta_1(:,1),age,lyr,inh,responsive_stat,'drift_spont Stationary');


%%%OS
tunedOSI_stat = responsive_stat & OSI (:,1)>=0.5;
%tunedOSI_run = responsive_run & OSI (:,2)>=0.5; % firing rate (responsiveness) criteria for whether cells enter subsequent statistical analysis
%tunedOSI_either = responsive_either & OSI (:,1)>=0.5;

layerAgeCDF(OSI(:,1),age,lyr,inh,responsive_stat,'OSI all Stationary');
%layerAgeCDF(OSI(:,2),age,lyr,inh,tunedOSI_run,'OSI all running'); 

layerAgePlot(OSI(:,1),age,lyr,inh,responsive_stat,'OSI all Stationary');
% layerAgePlot(OSI(:,2),age,lyr,inh,responsive_run,'OSI all Running');

layerAgeScatterMedian(OSI(:,1),age,lyr,inh,responsive_stat,'OSI all Stationary');
% layerAgeScatterMedian(OSI(:,2),age,lyr,inh,responsive_run,'OSI all running');

% part 1: create orientations variable
ori = 0:30:330;
ori = ori * pi /180;




%%%polar plot figure
 use=age'==1 & tunedOSI_stat;
 
 mx = cellfun(@(x) max(x(:)),drift_Ori_Sf);
 max_all=max(mx);
 

drift_os=(drift_all(use,1).orientfreq_all);
a = arrayfun(@(x)(getfield(x,'thetatuning')),drift_all(use,1),'UniformOutput',false);

ma=cellfun(@max,a);

norm_data = {a} / ma;
figure
 
 for j=1:length(a); 
%    
%     
%     % plot the tuning function of the three neurons 
%    polar([ori ori(1)], [w(j,:) w(j,1)],'k')
     figure
    tuning = a{j};
    tuning(tuning<0 & abs(tuning)<0.5)=0;
    if sum(tuning<0)==0
        polar([ori ori(1)], [ a{j} a{j}(1)]);hold on
    end
 end 
%%%drift_os = squeeze(d.orient_tune(r,f,:));
       %%% h=polar((0:8)*pi/4,[drift_os' drift_os(1)],color{r});

        figure
        h=polar((0:8)*pi/4,[drift_os(1)' drift_os(1)],color{r});hold on

%%%OS tuning width
tuned_OS_w_stat = tunedOSI_stat & drift_theta_w(:,1) >=0.1 & drift_theta_w(:,1) <1.2;
%tuned_OS_w_run = tunedOSI_run & drift_theta_w(:,2) >=0.1 & drift_theta_w(:,2) <1.04; % .2 radians =  ~12deg and 1.04radians = ~ 60deg 

%tuned_OS_w_either = tunedOSI_either & drift_theta_w(:,1) >=0.1 & drift_theta_w(:,1) <1.2;

OS_w_stat = responsive_stat & drift_theta_w(:,1) >=0.1 & drift_theta_w(:,1) <1.04;
%OS_w_run = responsive_run & drift_theta_w(:,2) >=0.1 & drift_theta_w(:,2) <1.04;

layerAgeCDF(drift_theta_w(:,1),age,lyr,inh,tuned_OS_w_stat ,'OS tuning width Stationary');
% layerAgeCDF(drift_theta_w(:,2),age,lyr,inh,tuned_OS_w_run ,'OS tuning width Running');

layerAgePlot(drift_theta_w(:,1),age,lyr,inh,tuned_OS_w_stat,'OS tuning width Stationary');
% layerAgePlot(drift_theta_w(:,2),age,lyr,inh,tuned_OS_w_run,'OS tuning width Running');

layerAgeScatterMedian(drift_theta_w(:,1),age,lyr,inh,tuned_OS_w_either,'OS tuning width Stationary');
% layerAgeScatterMedian(drift_theta_w(:,2),age,lyr,inh,tuned_OS_w_run,'OS tuning width Running');


%%%DS
DSI_stat = responsive_stat & driftdsi(:,1)>=0 & driftdsi(:,1)<1.05;
%DSI_run =  responsive_run & driftdsi(:,2)>=0.05 & driftdsi(:,2)<1.05;
DSI_either = responsive_either & driftdsi(:,1)>=0 & driftdsi(:,1)<1.05;

top30prct_stat = prctile(driftdsi(:,1),70);
DSI_stat_top30 = responsive_stat & driftdsi(:,1)>= top30prct_stat & driftdsi(:,1)<1.05;

% top30prct_run = prctile(driftdsi(:,2),70);
% DSI_run_top30 = responsive_run & driftdsi(:,2)>= top30prct_run & driftdsi(:,2)<=1;

layerAgeCDF(driftdsi(:,1),age,lyr,inh,DSI_stat_top30,'DSI Stationary ');
%layerAgeCDF(driftdsi(:,2),age,lyr,inh,DSI_run_top30,'DSI Running');

layerAgeCDF(driftdsi(:,1),age,lyr,inh,DSI_stat,'tuned DSI Stationary ');
% layerAgeCDF(driftdsi(:,2),age,lyr,inh,DSI_run,'tuned DSI run ');

layerAgePlot(driftdsi(:),age,lyr,inh,DSI_stat_top30,'DSI Stationary ');
%layerAgePlot(driftdsi(:,2),age,lyr,inh,DSI_run_top30,'DSI Running');
 
layerAgePlot(driftdsi(:,1),age,lyr,inh,DSI_stat,'DSI Stationary ');
% layerAgePlot(driftdsi(:,2),age,lyr,inh,DSI_run,'DSI Running');

layerAgeScatterMedian(driftdsi(:,1),age,lyr,inh,DSI_stat_top30,'DSI Stationary');
% layerAgeScatterMedian(driftdsi(:,2),age,lyr,inh,tuned_DSI_run,'DSI Running');
% 
layerAgeScatterMedian(driftdsi(:,1),age,lyr,inh,DSI_stat,'DSI Stationary');
%layerAgeScatterMedian(driftdsi(:,2),age,lyr,inh,DSI_run_top30,'DSI Running');

%%%receptive field size(width in radians)
rfw_tuned_stat = tunedOSI_stat & rfw(:,1)<=16& rfw(:,1)>=1;
%rfw_tuned_run = tunedOSI_run & rfw(:,2)<=16 & rfw(:,2)>=1;

rfw_all_stat =  responsive_stat & rfw(:,1)<=16 & rfw(:,1)>=1;
%rfw_all_run =  responsive_run & rfw(:,2)<=16 & rfw(:,2)>=1;
rfw_all = responsive_either & rfw(:,1)<=16 & rfw(:,1)>=1;

layerAgeCDF(rfw(:,1),age,lyr,inh,rfw_all_stat,'RFsize Stationary ');
% layerAgeCDF(rfw(:,2),age,lyr,inh,rfw_all_limit_run,'RFsize Running ');

layerAgeCDF(rfw(:,1),age,lyr,inh,rfw_tuned_stat,'RFsize Stationary ');
 
% layerAgeCDF(rfw(:,1),age,lyr,inh,rfw_tuned_limit_stat,'RFsize Stationary');
% layerAgeCDF(rfw(:,2),age,lyr,inh,rfw_tuned_limit_run,'RFsize Running ');

layerAgePlot(rfw(:,1),age,lyr,inh,rfw_tuned_stat,'RFsize Stationary');
% layerAgePlot(rfw(:,2),age,lyr,inh,rfw_all_limit_run,'RFsize Running ');

layerAgePlot(rfw(:,1),age,lyr,inh,rfw_all_stat,'RFsize Stationary ');
% layerAgePlot(rfw(:,1),age,lyr,inh,rfw_tuned_limit_stat,'RFsize Stationary');
% layerAgePlot(rfw(:,2),age,lyr,inh,rfw_tuned_limit_run,'RFsize Running ');

layerAgeScatterMedian(rfw(:,1),age,lyr,inh,rfw_all_stat,'RFsize Stationary');
% layerAgeScatterMedian(rfw(:,2),age,lyr,inh,rfw_all_limit_run,'RFsize Running');
% 
layerAgeScatterMedian(rfw(:,1),age,lyr,inh,rfw_tuned_stat,'RFsize Stationary ');

% layerAgeScatterMedian(rfw(:,1),age,lyr,inh,rfw_tuned_limit_stat,'RFsize Stationary');
% layerAgeScatterMedian(rfw(:,2),age,lyr,inh,rfw_tuned_limit_run,'RFsize Running ');


%%%spatial frequency preference and bandwidth
 

driftwpref(driftwpref==0) = 0.005; %tranform SF pref = 0 to 0.005 in order to plot on log scale


 SF_pref_stat = responsive_stat & driftwpref(:,1) >=0.005 & driftwpref(:,1)<=0.40;
% SF_pref_run = responsive_run & driftwpref(:,2) >=0.005 & driftwpref(:,2)<=0.32;

% SF_pref_stat_SF = responsive_stat & driftwpref(:,1) >=0 & driftwpref(:,1)<=8;
% SF_pref_run_SF = responsive_run & driftwpref(:,2) >0 & driftwpref(:,2)<=.32;
% 
% SF_pref_stat = responsive_stat &  driftwpref(:,1)<=.32;
% SF_pref_run = responsive_run  & driftwpref(:,2)<=.32;
% 
% SF_pref_tuned_stat = tunedOSI_stat & driftwpref(:,1)<=.32;
% SF_pref_tuned_run = tunedOSI_run & driftwpref(:,2) <=.32;
% 
% SF_pref_tuned_stat_SF = tunedOSI_stat & driftwpref(:,1)>0 & driftwpref(:,1)<=.32;
% SF_pref_tuned_run_SF = tunedOSI_run & driftwpref(:,2)>0 & driftwpref(:,2) <=.32;

%SF_pref_stat = responsive_either & driftwpref(:,1)>=0.005;
%SF_pref_stat_either = responsive_either & driftwpref(:,1)>=1 & driftwpref(:,1)<=8;
% SF_pref_tuned_run = tunedOSI_run &  driftwpref(:,2) >=0.005 & driftwpref(:,2)<=0.32;

layerAgeCDF(driftwpref(:,1),age,lyr,inh,SF_pref_stat,'wpref Stationary');
%layerAgeCDF(driftwpref(:,2),age,lyr,inh,SF_pref_run,'wpref run');

%layerAgeCDF(driftwpref(:,1),age,lyr,inh,SF_pref_stat_either ,'wpref Stationary');
%layerAgeCDF(driftwpref(:,2),age,lyr,inh,SF_pref_tuned_run_SF,'wpref Running');

layerAgePlot(driftwpref(:,1),age,lyr,inh,SF_pref_stat,'wpref_stat'); 
%layerAgePlot(driftwpref(:,2),age,lyr,inh,SF_pref_run_SF,'wpref_run');

%layerAgePlot(driftwpref(:,1),age,lyr,inh,SF_pref_stat_either,'wpref_stat'); 
%layerAgePlot(driftwpref(:,2),age,lyr,inh,SF_pref_tuned_run,'wpref_run');

layerAgeScatterMedian(driftwpref(:,1),age,lyr,inh,SF_pref_stat_SF,'wpref_stat');
% layerAgeScatterMedian(driftwpref(:,2),age,lyr,inh,SF_pref_run,'wpref_run');


%%%SF_pref_bandwidth

driftwbw_1=size(driftwbw);
driftwbw_1=driftwbw;
driftwbw_1(logical(imag(driftwbw_1)))=4;


driftwbw_1(driftwbw_1<=0.01)=0.05;
driftwbw_1(driftwbw_1==4)=4.5;
driftwbw_1(driftwbw_1==0.05)=4;

% driftwbw_2 =size(driftwbw_1)
% driftwbw_2=log(driftwbw_1);

bw_stat =  responsive_stat %% & driftwbw_1(:,1)<=3.8; %% & driftwbw_2(:,1)<=1.5;


layerAgeCDF(driftwbw_1(:,1),age,lyr,inh,bw_stat,'SF_pref tuning width Stationary');
% layerAgeCDF(drift_theta_w(:,2),age,lyr,inh,tuned_OS_w_run ,'OS tuning width Running');

layerAgePlot(driftwbw_1(:,1),age,lyr,inh,bw_stat,'SF_pref tuning width Stationary');
% layerAgePlot(drift_theta_w(:,2),age,lyr,inh,tuned_OS_w_run,'OS tuning width Running');

layerAgeScatterMedian(driftwbw_1(:,1),age,lyr,inh,bw_stat,'SF_pref tuning width Stationary');
% layerAgeScatterMedian(drift_theta_w(:,2),age,lyr,inh,tuned_OS_w_run,'OS tuning width Running');


%%%simple cells F1F0
tunedOSI_F1F0_stat = tunedOSI_stat & driftF1F0(:,1) <2.2 ;
F1F0_stat = responsive_stat& driftF1F0(:,1) <2.2;

% tunedOSI_F1F0_run = tunedOSI_run & driftF1F0(:,2) <=2.1 ;
 %F1F0_run = responsive_run & driftF1F0(:,2) <=2.1;

 
layerAgeCDF(driftF1F0(:,1),age,lyr,inh,F1F0_stat,'F1F0 Stationary');% F1F ratio of all units regardless of selectivity
%layerAgeCDF(driftF1F0(:,2),age,lyr,inh,F1F0_run,'F1F0 Running');

layerAgeCDF(driftF1F0(:,1),age,lyr,inh,tunedOSI_F1F0_stat,'F1F0 Stationary');
% layerAgeCDF(driftF1F0(:,2),age,lyr,inh,tunedOSI_F1F0_run,'F1F0 Running');
 
layerAgePlot(driftF1F0(:,1),age,lyr,inh, F1F0_stat,'F1F0 Stationary');
%layerAgePlot(driftF1F0(:,2),age,lyr,inh, F1F0_run,'F1F0 Running');

layerAgePlot_frac_simple(driftF1F0(:,1),age,lyr,inh, F1F0_stat,'F1F0 Stationary');
%layerAgePlot_frac_simple(driftF1F0(:,2),age,lyr,inh, F1F0_run,'F1F0 Running');

layerAgePlot_frac_simple(driftF1F0(:,1),age,lyr,inh, tunedOSI_F1F0_stat,'F1F0 Stationary');
%layerAgePlot_frac_simple(driftF1F0(:,2),age,lyr,inh, F1F0_run,'F1F0 Running');

layerAgeScatterMedian(driftF1F0(:,1),age,lyr,inh, F1F0_stat,'F1F0 Stationary');
% layerAgeScatterMedian(driftF1F0(:,2),age,lyr,inh, F1F0_run,'F1F0 Running');

layerAgeScatterMedian(driftF1F0(:,1),age,lyr,inh, tunedOSI_F1F0_stat,'F1F0 Stationary');
% layerAgeScatterMedian(driftF1F0(:,2),age,lyr,inh, tunedOSI_F1F0_run,'F1F0 Running');


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



% <<<<<<< HEAD
% keyboard
layerAgePlot(wn_evoked'/30,age,lyr,inh,1,'wn evoked');
layerAgePlot(wn_spont'/30,age,lyr,inh,1,'wn spont');
layerAgeScatter(wn_spont',age,lyr,inh,1,'wn spont');
layerAgeScatter(wn_evoked',age,lyr,inh,1,'wn evoked');
% =======
%keyboard
layerAgePlot(wn_evoked/30,age,lyr,inh,1,'wn evoked');
layerAgePlot(wn_spont,age,lyr,inh,1,'wn spont');
layerAgeScatter(wn_spont,age,lyr,inh,1,'wn spont');
layerAgeScatter(wn_evoked,age,lyr,inh,1,'wn evoked');
% >>>>>>> 160c0cebef96e05de28c5d9f5854ed64fd9595dc

layerAgeScatter(peak(:,2),age,lyr,inh,1,'moving peak');
figure
hist(wn_evoked); xlabel('wn_evoked')
figure
%%% SF

figure
plot(wn_stop(:,1),wn_mv(:,1),'.'); axis equal; hold on
plot([0 20],[0 20],'g'); xlabel('stationary'); ylabel('moving')
title('spont')

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



figure
use = ~midnarrow;
plot(wn_stop(use,10)-wn_stop(use,1),wn_mv(use,10)-wn_mv(use,1),'b.'); axis equal;hold on
use = midnarrow;
plot(wn_stop(use,10)-wn_stop(use,1),wn_mv(use,10)-wn_mv(use,1),'r.'); axis equal;hold on
plot([0 20],[0 20],'k');xlabel('stationary'); ylabel('moving')
use = inh;
plot(wn_stop(use,10)-wn_stop(use,1),wn_mv(use,10)-wn_mv(use,1),'g.'); axis equal;hold on
title('evoked')

layerAgeActivity(wn_stop(:,10)-wn_stop(:,1),wn_mv(:,10)-wn_mv(:,1),age,lyr,inh,1,{'stationary','moving'},'evoked');


LFPmax = max(nanmax(moveLFP,[],3),[],2);
figure
hist(LFPmax)
moveLFP(:,:,59:61)=0.5*moveLFP(:,:,59:61);

LFPmax_Bar = max(nanmax(moveLFP_bar,[],3),[],2);
figure
hist(LFPmax_Bar)
moveLFP_bar(:,:,59:61)=0.5*moveLFP_bar(:,:,59:61);

layerAgePlotMv([wn_stop(:,1) wn_mv(:,1)],age,lyr,inh,1,'spont');

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


[ m e n] =layerAgePlotMv([wn_stop(:,10)-wn_stop(:,1),wn_mv(:,10)-wn_mv(:,1)],age,lyr,inh,1,'evoked');

[ m e n] =layerAgePlotMv([peak(:,1),peak(:,2)],age,lyr,inh,peak(:,1)>1 & peak(:,2)>1,'drift evoked');

[ m e n] =layerAgePlotMv([driftspont(:,1),driftspont(:,2)],age,lyr,inh,1,'drift spont');

%plot the LFP power spectrum for the whole population during WN stim
figure
imagesc(squeeze(moveLFP(:,1,:)))

figure
imagesc(squeeze(moveLFP(:,2,:)))

figure
imagesc(squeeze(moveLFP(LFPmax<1*10^4,2,:)))

figure
imagesc(squeeze(moveLFP(LFPmax<1*10^4,1,:)))

%%bar stim LFPs

% figure
% imagesc(squeeze(moveLFP_bar(:,1,:)))

figure
imagesc(squeeze(moveLFP_bar(LFPmax_Bar<1*10^4,1,:)))

figure
imagesc(squeeze(moveLFP_bar(LFPmax_Bar<1*10^4,2,:)))

ly = 5;
figure
imagesc(squeeze(moveLFP_bar(age'==1 & lyr <ly & LFPmax_Bar<1*10^4,1,:)))

figure
imagesc(squeeze(moveLFP_bar(age'==1 & lyr <ly & LFPmax_Bar<1*10^4,2,:)))

%%filter LFP signals
ly =5
H = fspecial('disk', 1);
sm_moveLFP = imfilter(squeeze(moveLFP),H);

figure
imagesc(squeeze(sm_moveLFP(age'==2 & lyr<= ly & LFPmax<2.5*10^4 ,2,:)))
hold on
% title 'unfilterd'

title ' Bar running average filter [1 10 ]'


figure
imagesc(squeeze(sm_moveLFP_bar(age'==2 & lyr <ly  ,2,:)))
%broken down by layer and group
ly =5;
figure
imagesc(squeeze(moveLFP_bar(age'==2 & lyr <ly & LFPmax_Bar<2*10^4 ,2,:)))



figure
imagesc(squeeze(moveLFP_bar(LFPmax_Bar<2*10^4,1,:)))



%plot the LFP power running vs. stationary for EO1 and adult populations 
figure
for ly = 2:5
    plot(1:120,squeeze(nanmedian(sm_moveLFP(age'==1 & lyr==ly & LFPmax<2*10^4,1,:),1)),'Linewidth',1+ly/3);
hold on
plot(1:120,squeeze(nanmedian(sm_moveLFP(age'==1& lyr==ly& LFPmax<2*10^4,2,:),1)),'g','Linewidth',1+ly/3);
end
title('median EO')

figure
for ly = 2:6
    plot(1:120,squeeze(nanmean(moveLFP(age'==1 & lyr==ly & LFPmax<2*10^4,1,:),1)),'Linewidth',1+ly/3);
hold on
plot(1:120,squeeze(nanmean(moveLFP(age'==1& lyr==ly& LFPmax<2*10^4,2,:),1)),'g','Linewidth',1+ly/3);
end
title('mean EO')


figure
for ly = 2:6
    plot(1:120,squeeze(nanmedian(moveLFP(age'==2 & lyr==ly & LFPmax<2*10^4,1,:),1)),'Linewidth',1+ly/3);
hold on
plot(1:120,squeeze(nanmedian(moveLFP(age'==2 & lyr==ly & LFPmax<2*10^4,2,:),1)),'g','Linewidth',1+ly/3);
end
title('median adult')

figure
for ly = 2:6
    plot(1:120,squeeze(nanmean(moveLFP(age'==2 & lyr==ly & LFPmax<2*10^4,1,:),1)),'Linewidth',1+ly/3);
hold on
plot(1:120,squeeze(nanmean(moveLFP(age'==2 & lyr==ly & LFPmax<2*10^4,2,:),1)),'g','Linewidth',1+ly/3);
end
title('mean adult')

%%%LFP bars pop avg
%plot the LFP power running vs. stationary for EO1 and adult populations 
figure
for ly = 2:6
    plot(1:120,squeeze(nanmedian(moveLFP_bar(age'==1 & lyr==ly & LFPmax_Bar<2*10^4,1,:),1)),'Linewidth',1+ly/3);
hold on
plot(1:120,squeeze(nanmedian(moveLFP_bar(age'==1& lyr==ly& LFPmax_Bar<2*10^4,2,:),1)),'g','Linewidth',1+ly/3);
end
title('median EO')

figure
for ly = 2:6
    plot(1:120,squeeze(nanmean(moveLFP_bar(age'==1 & lyr==ly & LFPmax_Bar<2*10^4,1,:),1)),'Linewidth',1+ly/3);
hold on
plot(1:120,squeeze(nanmean(moveLFP_bar(age'==1& lyr==ly& LFPmax_Bar<2*10^4,2,:),1)),'g','Linewidth',1+ly/3);
end
title('mean EO')


figure
for ly = 2:6
    plot(1:120,squeeze(nanmedian(moveLFP_bar(age'==2 & lyr==ly & LFPmax_Bar<2*10^4,1,:),1)),'Linewidth',1+ly/3);
hold on
plot(1:120,squeeze(nanmedian(moveLFP_bar(age'==2 & lyr==ly & LFPmax_Bar<2*10^4,2,:),1)),'g','Linewidth',1+ly/3);
end
title('median adult')

figure
for ly = 2:6
    plot(1:120,squeeze(nanmean(moveLFP_bar(age'==2 & lyr==ly & LFPmax_Bar<2*10^4,1,:),1)),'Linewidth',1+ly/3);
hold on
plot(1:120,squeeze(nanmean(moveLFP_bar(age'==2 & lyr==ly & LFPmax_Bar<2*10^4,2,:),1)),'g','Linewidth',1+ly/3);
end
title('mean adult')
