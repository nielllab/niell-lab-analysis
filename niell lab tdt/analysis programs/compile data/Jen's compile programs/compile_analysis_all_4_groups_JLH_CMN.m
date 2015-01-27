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
        
        afiles = { 'Good recordings\Adults\9_12_13\rec1\analysis_9_12_13_adult_rec1.mat',...
             'Good recordings\Adults\11_15_13\rec1\analysis_11_15_13_adult_rec1.mat',...
             'Good recordings\Adults\1month_2month old\7_18_13_adult\analysis_07_18_13_cluster_rec1.mat',...
            'Good recordings\Adults\5_25_13_Adultmale3mo\Adult B\analysis_5_25_13_B_rec1_adult.mat',...
            'Good recordings\Adults\1month_2month old\4_22_13\analysis_adult.mat',...
            'Good recordings\Adults\7_19_13_adult\analysis_cluster_adult_rec2_7_19_13.mat',...
            'Good recordings\Adults\7_25_13_deep\MouseB\Rec1\analysis_7_25_13_mouseB_rec1.mat',...
             'Good recordings\Adults\11_13_13\rec1\analysis_adult_11_13_13_rec1.mat',...
              'Good recordings\Adults\11_14_13\rec1\analysis_11_14_13_adult_rec1.mat',...
            'Good recordings\Adults\11_14_13\rec2\analysis_11_14_13_adult_rec2.mat'};
%             'Good recordings\Adults\9_12_13\rec2\analysis_9_12_13_adult_rec2.mat',...
%             'Good recordings\Adults\11_11_13\rec1\analysis_11_11_13_adult_rec1.mat',... %%% no wn
%             'Good recordings\Adults\11_11_13\rec2\analysis_11_11_13_adult_rec2.mat',...
             
%             'Good recordings\Adults\11_15_13\rec2\analysis_11_15_13_adult_rec2.mat',... %% no wn
 %           'Good recordings\Adults\1month_2month old\9_25_13\rec1\analysis_9_25_13_P31_rec1.mat',...
            
%             'Good recordings\Adults\7_25_13_deep\MouseB\Rec2\analysis_cluster_data_07_25_13_mouseB_adult_rec2.mat',...
%             'Good recordings\Adults\1month_2month old\9_25_13\rec1\analysis_9_25_13_P31_rec1.mat',...
       %      'Good recordings\Adults\11_13_13\rec1\analysis_adult_11_13_13_rec1.mat',...
  %          'Good recordings\Adults\11_13_13\rec2\analysis_11_13_13_rec2.mat',...
          
            
%     elseif dataset==2
%         
%         afiles = {'Good recordings\EO7_EO9\9_6_13_EO11\analysis_9_6_13_EO11_rec1.mat',...
%               'Good recordings\EO7_EO9\11_27_13_EO7\rec1\analysis_11_27_13_EO7_rec1.mat',...
%             'Good recordings\EO7_EO9\P24_9_11_14\analysis_9_11_14_rec1.mat',...
%              'Good recordings\EO7_EO9\9_30_14_EO9\analysis_9_30_14_EO9_rec1.mat',...
%              'Good recordings\EO7_EO9\8_6_13_EO7\analysis_EO7_rec1_8_6_13.mat',...
%               'Good recordings\EO7_EO9\8_6_13_EO7\rec2\analysis_8_6_13_rec2_EO7.mat',...
%               'Good recordings\EO7_EO9\9_2_13_EO7\rec2\analysis_9_2_13_EO9_rec2.mat',...
%               'Good recordings\EO7_EO9\9_6_13_EO11\rec2\analysis_9_6_13_EO11_rec2.mat',...
%               'Good recordings\EO7_EO9\11_27_13_EO7\rec2\analysis_11_27_13_EO7_rec2.mat',...
%               'Good recordings\EO7_EO9\05_10_13_EO7\Record1_upper\analysis_1.mat',...
%               'Good recordings\EO7_EO9\05_10_13_EO7\Record_1_deeper\analysis_EO8_deeper.mat',...
%               'Good recordings\EO7_EO9\P26_9_15_14\analysis_9_15_14.mat',...
%              'Good recordings\EO7_EO9\P24_9_11_14\rec2\analysis_9_11_14_rec2.mat',...
%             'Good recordings\EO7_EO9\10_2_14_EO5\rec1\analysis_10_2_14_rec1_EO5.mat'};
%     elseif dataset ==3
%        
%         afiles = { 'Good recordings\EO3_EO4\7_5_13P16_EO4\mouseB\analysis_07_5_13_cluster_7_5_13_mouseB_945uM_analysis.mat',...
%              'Good recordings\EO3_EO4\8_2_13_EO3\analysis_EO3_8_2_13_rec1.mat',...
%              'Good recordings\EO3_EO4\9_10_13_EO3\rec1\analysis_9_10_13_EO3_rec1.mat',...
%              'Good recordings\EO3_EO4\9_10_13_EO3\rec2\analysis_9_10_13_EO3_rec2.mat',...
%                          'Good recordings\EO3_EO4\EO4_9_25_14_rec1\analysis_9_25_14_EO4_rec1.mat',...
%                          'Good recordings\EO3_EO4\5_6_13_EO3\mouseC\analysis_mouseC_EO3_5_6_13.mat',...
%              'Good recordings\EO3_EO4\11_25_13_EO4\rec1\analysis_11_25_13_EO4_rec1.mat',...
%              'Good recordings\EO3_EO4\11_25_13_EO4\rec2\analysis_11_25_13_EO4_rec2.mat',...
%              'Good recordings\EO3_EO4\10_1_14_EO4\analysis_10_1_14_rec1_EO4_Altered.mat',...
%             'Good recordings\EO3_EO4\10_1_14_EO4\rec2\analysis_10_1_14_rec2_EO4.mat',...
%             'Good recordings\EO3_EO4\EO3_9_24_14\analysis_9_24_14_EO3.mat',...
%             'Good recordings\EO3_EO4\EO4_9_25_14_rec1\rec2\analysis_9_25_14_EO4_rec2.mat'}; %% no wn
          
     
    elseif dataset ==2
        afiles = {'Good recordings\EO1_EO2\8_7_13_EO1\rec1_full_clustering\analysis_8_7_13_EO1_rec1_more_strigent.mat',...
            'Good recordings\EO1_EO2\8_8_13_EO2\rec1_full_clustering\analysis.mat',...
            'Good recordings\EO1_EO2\5_22_13_EO1\analysis_rec1_A_5_22_13_strict_selection.mat',...
            'Good recordings\EO1_EO2\7_17_13_EO1\Analysis_7_17_13_cluster_7_17_13_EO1.mat',...
            'Good recordings\EO1_EO2\7_17_13_EO1\Rec2\Analysis_7_17_13_rec2_EO1.mat',...
            'Good recordings\EO1_EO2\9_9_13_EO1\rec2\analysis_9_9_13_EO1_rec2.mat',...
            'Good recordings\EO1_EO2\10_1_13_EO2\rec1\analysis_10_1_13_EO2_rec1.mat',...
            'Good recordings\EO1_EO2\11_21_13_EO1_mouseA\analysis_11_21_13_mouseA_EO1.mat',...
            'Good recordings\EO1_EO2\11_23_13_EO2\analysis_11_23_13_EO2.mat'};
      end  
           % 'Good recordings\EO1_EO2\9_9_13_EO1\rec1\analysis_9_9_13_EO1_rec1.mat',...
            %'Good recordings\EO1_EO2\9_30_13_EO1\rec1\analysis_9_30_13_EO1_rec1.mat',...
            %'Good recordings\EO1_EO2\9_30_13_EO1\rec2\analysis_9_30_13_EO1_rec2.mat',...
            
            %'Good recordings\EO1_EO2\11_20_13_EO0\rec2\analysis_11_20_13_EO0_rec2.mat',...
           
           % 'Good recordings\EO1_EO2\11_26_13_mouseB_EO1\analysis_11_26_13_mouseB_EO1.mat',...
           % 'Good recordings\EO1_EO2\11_27_13_mouseB_EO2\analysis_11_27_13_EO2_mouseB.mat'}; %%% tg
      %           'Good recordings\EO1_EO2\11_20_13_EO0\rec1\analysis_11_20_13_rec1.mat',...
%             'Good recordings\EO1_EO2\11_26_13_mouseA_EO1\rec1\analysis_11_26_13_EO1_mouseA_rec1.mat',...
%             'Good recordings\EO1_EO2\11_26_13_mouseA_EO1\rec2\analysis_11_26_13_EO1_rec2.mat',...
%             'Good recordings\EO1_EO2\8_7_13_EO1\rec2_full_clustering\analysis_8_7_13_rec2.mat',...
 %       'Good recordings\EO1_EO2\4_29_13_EO1\mouseC\analysis_4_29_13_C.mat',...
   
    
    
    for i = 1:length(afiles)
   
        clear params
        clear wn wn_movement
        clear LFP_movement
        clear bars wave_all
        clear rf_width
        clear locomotion
        load([apath afiles{i}]);
%         
%        clusterfilename
%        afiles{i}
%     
%         if exist(clusterfilename,'file')
%             clusterFile = clusterfilename;
%         elseif exist(clusterfilename((length(pname)+1):end),'file')
%             clusterFile = clusterfilename((length(pname)+1):end);
%         elseif exist([clusterfilename((length(pname)+1):end) '.mat'],'file')
%             clusterFile = [clusterfilename((length(pname)+1):end) '.mat'];
%         else
%             [fname pname] = uigetfile('*.mat','cluster file');
%             clusterFile = fullfile(pname,fname);
%             clusterfilename = clusterFile;
%             if fname~=0
%                 save([apath afiles{i}],'clusterfilename','-append');
%             end
%         end
%         clusterFile
%         try
%             load(clusterFile,'wave_all');
%         catch
%             display('no cluster file')
%         end
%         
%     
         
        
        
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
        
        %size(rf_width)
        
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
        
        driftA1_th(cellrange,:)= field2array(drift_th,'A1');
        driftA2_th(cellrange,:)=field2array(drift_th,'A2');
        driftB_th(cellrange,:)= field2array(drift_th,'B');
        drift_theta_w_th(cellrange,:)=field2array(drift_th,'thetawidth');
        drift_theta_th(cellrange,:)=field2array(drift_th,'theta');
        
        
        
        driftspont(cellrange,:) = field2array(drift,'spont');
        driftwpref(cellrange,:) = field2array(drift,'wpref');
        driftwbw(cellrange,:) = field2array(drift,'bw') ;
        
%        driftF1F0(cellrange,:) = field2array(drift,'F1')./field2array(drift,'F0');
        driftF0(cellrange,:) = field2array(drift,'F0');
        %       driftorientfreq_all(cellrange,:)=field2array(drift, 'orientfreq_all');

        driftlayer =  field2array(drift,'layer');
        lyr(cellrange,:) = driftlayer(:,1);
        cvDSI(cellrange,:) = field2array(drift,'dsi'); %%also circular variance measure change to "cv_dsi and cv_osi" in new compile programs
        cvOSI(cellrange,:)=field2array(drift,'osi');
        
        SNR(cellrange,:)=field2array(drift_sig_noise,'signoise');
        SNR_err(cellrange,:)=field2array(drift_sig_noise,'signoise_SE');
        %driftOri(cellrange,:) = field2array(drift,'orientfreq_all');
        
        if exist('bars');

        bar_spont(cellrange,:)=field2array(bars,'spont');
       
        else
     %   bar_spont(cellrange,:)= NaN;
        
        end
        
       
        
%         clear meanwaves snr stdwaves
%         
%         for c=1:length(cells);
%             %%% get SNR
%             tet =ceil(cells(c,1)/4);
%             
%             if exist('wave_all','var')
%                 wvall = wave_all{tet};
%                 wvclust = wvall(find(idx_all{(tet-1)*4+1}==cells(c,2)),:,:);
%                 
%                 amps =squeeze(min(wvclust(:,5:10,:),[],2));
%                 mn = abs(nanmean(amps));
%                 stdev = nanstd(amps);
%                 [y ind] = max(mn);
%                 snr(c) = mn(ind)/stdev(ind);
%                 
%                 meanwaves(c,:,:) = squeeze(nanmean(wvclust,1));
%                 stdwaves(c,:,:) = squeeze(nanstd(wvclust,[],1));
%             else
%                 meanwaves=NaN;
%                 snr=NaN
%                 stdwaves=NaN
%             end
%             
%             
%         end
        
%         SNRall(cellrange)=snr;
%         meanWavesAll(cellrange,:,:) = meanwaves;
%         stdWavesAll(cellrange,:,:) = stdwaves;
        
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
% %         STA_nx(cellrange)= NaN;
% %         STA_ny(cellrange)=  NaN;
% %         STA_phase(cellrange)= NaN;
% %         STA_sigx(cellrange)= NaN;
% %         STA_sigy(cellrange)= NaN; 
% %         STA_exp_var(cellrange)=NaN;
% 
%         end
        

        
        %size(wvform)
        wvform(cellrange,:) = wv';
     
        %get firing rate at all measured orients and SF, put into an array:
        %12 rows(orientations) by 7 columns(SpatialFreqs) for each cell
        
        drift_Ori_Sf(cellrange,:) = arrayfun(@(x)(getfield(x,'orientfreq_all')),drift,'UniformOutput',false);
       % drift_all(cellrange,:)=drift';
      
       
       if exist('locomotion');
           Vel(i,dataset)=arrayfun(@(x)(getfield(x,'mouseV')),locomotion,'UniformOutput',false);
           
           
       end
       
       end %%% loop over adult vs EO
end

% figure
% plot(squeeze(meanWavesAll(1,:,1)))
% 
% 
% waveforms=figure
% x=[1:19]
% 
% SNR_1=age==1 & SNRall>8 & SNRall<6
% [x idx]=max(SNRall(SNR_1))
% 
% SNR_7 =SNRall(age==3);
% SNR_3=SNRall(age==2);
% SNR_1=SNRall(age==1);

% for i=1:4
%     clear m s
%     
% m=squeeze(meanWavesAll(25,:,i))
% s =squeeze(stdWavesAll(25,:,i))   
% figure(waveforms)
% subplot(4,1,i); hold on
% shadedErrorBar(x,m,s)
% end
%plot speed distributions

Aspeed=zeros;






figure
colorlist='bgrm'
for dev= 1:2
speedlist=[]; 
for i = 1:size(Vel,1);
      frac_loc(i,dev)=(sum(Vel{i,dev}<1))/length(Vel{i,dev});
    %if frac_loc(i,age)<0.82
       speedlist=[speedlist Vel{i,dev}];
    %end 
end
speedlist = speedlist*9.8/7;
speedlist(speedlist<0.25)=0.25;
x=nanmedian(speedlist>1)

speedlist = log10(speedlist);

%[f3 h3]=hist(speedlist,0.5:1:20);
[f3 h3]=hist(speedlist,-0.5:0.25:2);
f3=f3./sum(f3);

plot(h3,f3,'color',colorlist(dev));hold on
xlim([-0.5 2]); ylim([0 0.55])
set(gca,'Xtick',[-0.5 :0.5 : 1.5])
set(gca,'Xticklabel',{'0','1','3','10','30'})
xlabel('speed cm/sec')

end

legend('E01','Adult')
plot([log10(0.7) log10(0.7)],[0 0.55],'k--')








% ind = find(isnan(driftF1F0));
% driftF1F0(ind)=0;


%define excitatory cell types versus inhibotory types based on trough
%to peak versus troughdepth/peak height

EI = [size_speed_1];
opts = statset('display','final');
[idx,ctrs] = kmeans(EI,2,'distance','city','Replicates',5,'options',opts);
if sum(idx==1)<sum(idx==2)
    run = (idx==1);
else
   run = (idx==2);
end



figure
plot(alldata(age==1,5),alldata(age==1,6),'bo'); hold on
plot(alldata(age==2,5),alldata(age==2,6),'mo'); hold on
plot(alldata(age==3,5),alldata(age==3,6),'ro'); hold on
plot(alldata(age==4,5),alldata(age==4,6),'go'); hold on

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
%plot(alldata(age'==4,5),alldata(age'==4,6),'go');hold on
%plot(alldata(age'==3,5),alldata(age'==3,6),'ro');hold on
plot(alldata(age'==2,5),alldata(age'==2,6),'color',[.5,.5,.5 ]);hold on
plot(alldata(age'==1,5),alldata(age'==1,6),'bo');
title 'EO7 vs. adults'

figure
plot(alldata(age'==4,5),alldata(age'==4,6),'go');hold on
plot(alldata(age'==2,5),alldata(age'==2,6),'mo')
title 'EO4 vs. adults'

figure
plot(wvform(age'==4& ~inh,:)','color','k');hold on
plot(wvform(age'==4& ~inh,:)','color','k');hold on
plot(wvform(age'==2 & ~inh,:)','color','r');hold on
plot(wvform(age'==4& inh,:)','color','r');hold on
 
figure
plot(wvform(age'==2,:)','color',[.5 .5 .5]);hold on

figure
 y=wvform(age'==1 ,:);
 %y=y';
 %plot(1:19,y(2,:));
 
 x=1:19;
 

for i=1:length(y); 
  
    y1= y(i,:);
    
         p  = patchline(x,y1,'edgecolor','b','linewidth',2,'edgealpha',0.2);
end
hold on

% figure
% plot(wvform(age'==4 & inh,:)','color','g','linewidth',2);hold on
% plot(wvform(age'==3 & inh,:)','color','r','linewidth',2);hold on
% plot(wvform(age'==2 & inh,:)','color','c','linewidth',2);hold on
% plot(wvform(age'==1 & inh,:)','color','b','linewidth',2);hold on
% 
% figure
% plot(wvform(age'==4,:)','color','g','linewidth',2);hold on
% plot(wvform(age'==3,:)','color','r','linewidth',2);hold on
% plot(wvform(age'==2,:)','color','c','linewidth',2);hold on
% plot(wvform(age'==1,:)','color','b','linewidth',2);hold on


figure
plot(wvform(age'==1,:)','g');hold on
plot(wvform(age'==2,:)','m');
title 'EO3_4 vs. adults'



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
        [OSI(i,j) DSI(i,j) width(i,j) peak(i,j)] = calculate_tuning(driftA1(i,j),driftA2(i,j),driftB(i,j),drift_theta_w(i,j));
        
        
    end
end

for i = 1:size(driftA1_th,1)
    for j=1:size(driftA1_th,2)
        driftA1_th(i,j);
        driftA2_th(i,j);
        driftB_th(i,j);
        drift_theta_w_th(i,j);
        [OSI_th(i,j) DSI_th(i,j) width_th(i,j) peak_th(i,j)] = calculate_tuning(driftA1_th(i,j),driftA2_th(i,j),driftB_th(i,j),drift_theta_w_th(i,j));
        
        
    end
end


%variables created to sift through data conditionally

%%firing rate
responsive_run = peak(:,2)>=2; 
responsive_stat = peak(:,1)>=2 & peak(:,2)>=2;  % firing rate (responsiveness) criteria for whether cells enter subsequent statistical analysis
responsive_stat_th = peak_th(:,1)>=1.8 & peak_th(:,2)>=2;  % firing rate (responsiveness) criteria for whether cells enter subsequent statistical analysis

responsive_SNR=SNR_err(:,1)>=2 & SNRall(:,1)>=2 
good_units=SNRall(:,1)>=2;

d= ~inh & age'==2 & lyr<=3 & responsive_stat;
g= ~inh & age'==1 & lyr<=3 & responsive_stat;
e=age'==1 & inh & responsive_stat;
f=age'==1 & inh & responsive_stat;
sum(d)
sum(g)
sum(e)

[p h]=ranksum(peak(d,1),peak(g,1))
% peak1=peak;
% peak1(peak1 < 0)=0.01; %%transforms all negative to 0.01 and makes that an equivelent class

evoked_stat= peak(:,1)>=2.0  & peak(:,1)<=80;  %OSI(:,1)>0.5 
evoked_run= responsive_run & peak(:,2)<=80;

all=peak(:,1)>=2;
%layerAgeCDF(peak1(:,1),age,lyr,inh,evoked_stat,'peak');
layerAgePlot(peak(:,1),age,lyr,inh,evoked_stat,'peak');
layerAgePlot_frac_responsive(peak(:,1),age,lyr,inh,responsive_SNR,good_units,'peak response');
layer_age_line(peak(:,1),age,lyr,inh,evoked_stat,'drift_peak FR Stationary')




figure

plot(peak(responsive_SNR & responsive_stat ,1),SNR(responsive_SNR& responsive_stat ,1),'mo')

i=nanmedian(peak1(inh & age'==4 & peak1(:,1)>=2))
%%%ratio of running peak to stat peak
layerAgePlot_ratio_jlh(peak_th(:,1),peak_th(:,2),age,lyr,inh,responsive_stat_th,{'run','stat'},'run vs stat');


%driftspont1=driftspont(:,1);
spont_stat_drift =  peak(:,1)>=2 & peak(:,1)<=80 & driftspont(:,1)>=0.045  ;
layerAgePlot_ratio_jlh(driftspont(:,1),driftspont(:,2),age,lyr,inh,spont_stat_drift ,{'run','stat'},'run vs stat');

dr1=driftspont(:,1)
dr2=driftspont(:,2)
[p h]=kstest2(dr1(age'==1& spont_stat_drift),dr2(age'==4&spont_stat_drift))

%%firing rate modulation

figure
plot(peak(find(~inh & responsive_both & age'==4),1),peak(find(~inh& responsive_both & age'==4),2),'go');hold on
plot([0 60],[0 60])
xlabel ('stationary firing rate')
ylabel ('locomotion firing rate')

%%by layers

figure
plot(peak(find(~inh & responsive_both & age'==1 & lyr<=3),1),peak(find(~inh & responsive_both & age'==1 & lyr<=3),2),'ko');hold on
hold on
plot(peak(find(~inh & responsive_both & age'==1 & lyr==4),1),peak(find(~inh & responsive_both & age'==1 & lyr==4),2),'ro');hold on
hold on
plot(peak(find(~inh & responsive_both & age'==1 & lyr==5),1),peak(find(~inh & responsive_both & age'==1 & lyr==5),2),'co');hold on
hold on
plot(peak(find(~inh & responsive_both & age'==1 & lyr==6),1),peak(find(~inh & responsive_both & age'==1 & lyr==6),2),'yo');hold on
hold on
plot([0 60],[0 60])
xlabel ('stationary firing rate')
ylabel ('locomotion firing rate')
legend ('layer2/3', 'layer4','layer5', 'layer6')
%[p h]=ranksum(exct,exct_run)

figure
plot(exct,exct_run,'go'); hold on
plot([0 60], [0 60])
plot(exct_E,exct_run_E,'bo'); hold on
title 'layer6'
%%conduct simple ranksum comparison of medians from two groups directly
 x1=responsive_stat & ~inh & age'==1 & lyr<=3;
 y1=responsive_stat & ~ inh & age'==4 & lyr<=3;

 [p h]=ranksum(peak1_stat(x1),peak1_stat(y1))



layerAgePlot(rfw(:,1),age,lyr,inh,responsive_stat,'RF size');
layer_age_line(rfw(:,1),age,lyr,inh,responsive_stat,'RF size');

peak1_stat=peak(:,1);
age_peak=age(peak(:,1)>=2);
lyr_peak=lyr(peak(:,1)>=2);
g={age_peak';lyr_peak};
p_peak=anovan(peak1_stat(peak(:,1)>=2),g, 'interaction');


 %%%%spontaneous rates during drift

driftspont1=driftspont(:,1);
spont_stat_drift =  peak(:,1)>=2 & driftspont(:,1)>=0.00  ;

%layerAgeCDF(driftspont1(:,1),age,lyr,inh,spont_stat_drift,'drift_spont Stationary');
layerAgePlot(driftspont(:,1),age,lyr,inh,spont_stat_drift,'drift_spont Stationary');
layer_age_line(driftspont(:,1),age,lyr,inh,spont_stat_drift,'drift_spont Stationary');

x1=responsive_stat & ~inh & age'==4 & driftspont(:,1)>=0 & lyr==5;
y1=responsive_stat & ~ inh & age'==4 & driftspont(:,1)>=0 &lyr==4;

 [p h]=ranksum(driftspont(x1,1),driftspont(y1,1))

%driftspont1(driftspont1<=.9)=1; %%%transforms zero and close to zero values to 1 for ease of taking log10
%driftspont1=log10(driftspont1);
bar_spont_1=bar_spont(:,1);
resp_stat_spont= responsive_stat & bar_spont_1(:,1)>=0.045;
%layerAgeScatterMedian(bar_spont_1(:,1),age,lyr,inh,resp_stat_spont,'bar_spont Stationary');
layerAgeCDF(bar_spont_1(:,1),age,lyr,inh,responsive_stat,'drift_spont Stationary');
layerAgePlot(bar_spont_1(:,1),age,lyr,inh,responsive_stat,'drift_spont Stationary');
p_peak=anovan(bar_spont_1(peak(:,1)>=2),g, 'interaction');

%ranksum(bar_spont_1(resp_stat_spont & age'==1 & lyr==4),bar_spont_1(resp_stat_spont & age'==2 & lyr==4))
%ranksum(driftspont1(spont_stat_drift & age'==1 & lyr==3),driftspont1(spont_stat_drift & age'==2 & lyr==3))

bar_spont_run=bar_spont(:,2);
resp_run_spont= responsive_run & bar_spont_run(:,1)>=0.045; 
%layerAgeScatterMedian(bar_spont_1(:,1),age,lyr,inh,resp_stat_spont,'bar_spont Stationary');
layerAgePlot(bar_spont_run,age,lyr,inh,resp_run_spont,'bar_spont runnig');
p_peak=anovan(bar_spont_run(peak(:,1)>=2),g, 'interaction');


%%%prefered orientation
drift_theta_1=(drift_theta*180)/pi;
drift_theta_1(drift_theta_1>330)=0;
%drift_theta_1(drift_theta_1>330)=0;


% clear top30prct_OSI OSI_stat_top30
% 
% top50prct_OSI = prctile(OSI(:,1),50);
% OSI_stat_top50 = responsive_either & OSI(:,1)>= top50prct_OSI; 
% 
% top50prct_OSI_run = prctile(OSI(:,2),50);
% OSI_run_top50 = responsive_run & OSI(:,2)>= top50prct_OSI_run;

responsive_stat=peak(:,1)>=2;
tunedOSI_stat = responsive_stat & OSI (:,1)>=0.5; %use=age'==2 & OSI (:,1)>=0.45 & lyr<=3;
tunedOSI_run = responsive_run & OSI (:,2)>=0.5; % firing rate (responsiveness) criteria for whether cells enter subsequent statistical analysis
tunedOSI_either = responsive_either & OSI (:,2)>=0.5|OSI (:,1)>=0.5 ;


%%preferred oprientation 
layerAgePlot_pref_Orient(drift_theta_1(:,1),age,lyr,inh,tunedOSI_stat ,{'pref Orient' 'Prct total'},'Prefered Orientation Stationary');

%layerAgePlot_frac_OS_DS(peak1_stat,age,lyr,inh,OS,c_DS,'OS that as c_DS'); 

%%%OS

layerAgeCDF(OSI(:,1),age,lyr,inh,responsive_stat,'OSI all Stationary');
%layerAgeCDF(OSI(:,2),age,lyr,inh,tunedOSI_run,'OSI all running'); 

layerAgePlot(OSI(:,1),age,lyr,inh,responsive_stat,'OSI all Stationary');
% layerAgePlot(OSI(:,2),age,lyr,inh,responsive_run,'OSI all Running');
layer_age_line(OSI(:,1),age,lyr,inh,responsive_stat,'OSI all Stationary');
%layerAgeScatterMedian(OSI(:,1),age,lyr,inh,responsive_stat,'OSI all Stationary');
% layerAgeScatterMedian(OSI(:,2),age,lyr,inh,responsive_run,'OSI all running');
o=OSI(:,1)
d=DSI(:,1)
clear x1 x2
x1=o(age'==4 & lyr==5 & responsive_stat);
x2=o(age'==4 & lyr==6 & responsive_stat);
[p h]= ranksum(x1,x2)

nanmedian(x1)
semedian(x1)
nanmedian(x2)
semedian(x2)
clear x1 x2
x1=d(age'==4 & lyr<=3 & responsive_stat);
x2=d(age'==3 & lyr<=3 & responsive_stat);
[p h]= ranksum(x1,x2)
y=sum(~isnan(x1))
y1=sum(~isnan(x2))
%%
% %%%DS
 DSI_stat =  responsive_stat & DSI(:,1)>=0 & DSI(:,1)<1.05;
 DSI_run =  DSI(:,2)>=0.5 & DSI(:,2)<1.05;
 DSI_either =  DSI(:,1)>=0 & DSI(:,1)<1.05;
% 
%  d=DSI(:,1)
% sf_H_ds=responsive_stat &  driftwpref(:,1) >=0.25 & driftwpref(:,1)<=0.40;
%  sf_L_ds=responsive_stat &  driftwpref(:,1) <0.25 & driftwpref(:,1)<=0.40;

%layerAgeCDF(cvDSI(:,1),age,lyr,inh,DSI_stat_top50,'DSI Stationary ');
%layerAgeCDF(cvDSI(:,1),age,lyr,inh,DSI_stat,'tuned DSI Stationary ');
 
layerAgePlot(DSI(:,1),age,lyr,inh,DSI_stat,'DSI Stationary ');
 
layer_age_line(cvDSI(:,1),age,lyr,inh,DSI_stat_top50,'DSI top 50%');
layer_age_line(cvDSI(:,1),age,lyr,inh,DSI_stat,'DSI ');




 OS=OSI(:,1)>=0.5 & responsive_stat;
 c_DS=cvDSI(:,1)>=0.5 & responsive_stat;
 OS_c_DS=OSI(:,1)>=0.5 & cvDSI(:,1)>=0.5 & responsive_stat;

 layerAgePlot_frac_OS_DS(peak(:,1),age,lyr,inh,OS,c_DS,'OS that as c_DS');

% part 1: create orientations variable
ori = 0:30:330;
ori = ori * pi /180;

%%%polar plot figure
 % use=age'==1 & responsive_stat & lyr==6;
 clear use
 %use=age'==2 & tunedOSI_stat  & DSI_stat_top50; %cvDSI(:,1)>=0.45
 use=age'==1 &  tunedOSI_stat & lyr<=3;% & DSI_stat_top50;
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
tuned_OS_w_stat = tunedOSI_stat & drift_theta_w(:,1) >=0.17 & drift_theta_w(:,1) <1.2;
%tuned_OS_w_run = tunedOSI_run & drift_theta_w(:,2) >=0.1 & drift_theta_w(:,2) <1.04; % .2 radians =  ~12deg and 1.04radians = ~ 60deg 
%tuned_OS_w_either = tunedOSI_either & drift_theta_w(:,1) >=0.1 & drift_theta_w(:,1) <1.2;

OS_w_stat=drift_theta_w(:,1) >=0.1 & drift_theta_w(:,1) <1.2;


OS_w_stat = responsive_stat & drift_theta_w(:,1) >=0.1 & drift_theta_w(:,1) <1.04;
%OS_w_run = responsive_run & drift_theta_w(:,2) >=0.1 & drift_theta_w(:,2) <1.04;

layerAgeCDF(drift_theta_w(:,1),age,lyr,inh,tuned_OS_w_stat ,'OS tuning width Stationary');
layerAgeCDF(drift_theta_w(:,1),age,lyr,inh,OS_w_stat ,'OS tuning width Stationary');

% layerAgeCDF(drift_theta_w(:,2),age,lyr,inh,tuned_OS_w_run ,'OS tuning width Running');

layerAgePlot(drift_theta_w(:,1),age,lyr,inh,OS_w_stat,'OS tuning width Stationary');
% layerAgePlot(drift_theta_w(:,2),age,lyr,inh,tuned_OS_w_run,'OS tuning width Running');

layer_age_line(drift_theta_w(:,1),age,lyr,inh,OS_w_stat,'OS tuning width');
% layerAgeScatterMedian(drift_theta_w(:,1),age,lyr,inh,tuned_OS_w_either,'OS tuning width Stationary');
% layerAgeScatterMedian(drift_theta_w(:,2),age,lyr,inh,tuned_OS_w_run,'OS tuning width Running');

drift=drift_theta_w(:,1);
 
x1=responsive_stat & drift_theta_w(:,1) >=0.1 & drift_theta_w(:,1) <1.04 & ~inh & age'==1 & lyr==6;
 y1=responsive_stat & drift_theta_w(:,1) >=0.1 & drift_theta_w(:,1) <1.04 & ~inh & age'==4 & lyr==6;
 [p h]=ranksum(drift(x1),drift(y1))

 
 %%spatial frequency preference and bandwidth
 
driftwpref(driftwpref==0) = 0.005; %tranform SF pref = 0 to 0.005 in order to plot on log scale


 SF_pref_stat = peak(:,1)>=2 & driftwpref(:,1) >=0.004 & driftwpref(:,1)<=0.40;
% SF_pref_run = responsive_run & driftwpref(:,2) >=0.005 & driftwpref(:,2)<=0.32;


layerAgeCDF(driftwpref(:,1),age,lyr,inh,SF_pref_stat,'wpref Stationary');
%layerAgeCDF(driftwpref(:,2),age,lyr,inh,SF_pref_run,'wpref run');

layerAgePlot(driftwpref(:,1),age,lyr,inh,SF_pref_stat,'wpref_stat'); 
layer_age_line(driftwpref(:,1),age,lyr,inh,SF_pref_stat,'SF_pref Stationary ');
layerAgePlot_SF_line(driftwpref(:,1),age,lyr,inh,SF_pref_stat,'SF_pref Stationary ');

all=peak(:,1);
responsive_stat=peak(:,1)>=2;
layerAgePlot_frac_responsive(driftwpref(:,1),age,lyr,inh,SF_pref_stat,'prct responsive');

driftwpref_stat=driftwpref(:,1)

age_peak=age(SF_pref_stat);
lyr_peak=lyr(SF_pref_stat);
g={age_peak';lyr_peak};
p_wpref=anovan(dift,g, 'interaction');
clear A

A={driftwpref_stat(SF_pref_stat & age'==4 & lyr<=6),driftwpref_stat(SF_pref_stat & age'==1&lyr<=6)};
[h p]= kstest2(A{1},A{2})
clear x1 y1
x1=SF_pref_stat & age'==1 & lyr==6
y1=SF_pref_stat & age'==4 & lyr==6

[p h]=ranksum(driftwpref_stat(x1),driftwpref_stat(y1))
d=sum(x1)
g=sum(y1)

figure
[f,x]=hist(driftwpref_stat(age'==1 & lyr<=6 & ~inh & SF_pref_stat),0:0.03:0.32);
H1=bar(x,f/sum(f),'b');
title 'EO1 nx L4'
hold on

[f1,x1]=hist(driftwpref_stat(age'==4 & lyr<=6 & ~inh& SF_pref_stat),0:0.03:0.32);
H2=bar(x,f1/sum(f1),'g');
ch=get(H2,'child');
set(ch,'facea',.5)
title 'Eo1_vs adult _ny'
hold on

%%%SF_pref_bandwidth

driftwbw_1=driftwbw(:,1);
driftwbw_1(logical(imag(driftwbw_1)))=-1;
driftwbw_1=1.5*(driftwbw_1);
%driftwbw_1(driftwbw_1==-1.5)=10;
driftwbw_1(driftwbw_1<1)=8;
% driftwbw_1(driftwbw_1<=7 & driftwbw_1>=6)=10;
 driftwbw_1(driftwbw_1<=0.4 & driftwbw_1>=-.5)=9;
driftwbw_1(driftwbw_1>8.5)=nan


figure
[f,x]=hist(driftwbw_1(age'==1 & ~inh & SF_pref_stat),0:0.5:9);
H1=bar(x,f/sum(f),'b');
title 'EO1 nx L4'
hold on

[f1,x1]=hist(driftwbw_1(age'==4 & lyr<=6 & ~inh & SF_pref_stat),0:0.5:9);
H2=bar(x,f1/sum(f1),'g');
ch=get(H2,'child');
set(ch,'facea',.5)
title 'Eo1_vs adult '
hold on



[f,x]=hist(driftwbw_1(age'==1 & SF_pref_stat ),0:0.5:9)
[f1,x1]=hist(driftwbw_1(age'==4 & SF_pref_stat ),0:0.5:9)
bw=[(f/sum(f));(f1/sum(f1))];
figure
bar(x,bw',2)
title 'EO1 vs afult bandwidth SF pref'

[f2,x2]=hist(driftwbw_1(age'==1 ),0:0.25:5);
[f3,x3]=hist(driftwbw_1(age'==3 ),0:0.25:5)

B={driftwbw_1(SF_pref_stat & age'==1 & lyr<=6),driftwbw_1(SF_pref_stat & age'==4 & lyr<=6)};

[h p]= kstest2(B{1},B{2})
[h p]= ranksum(B{1},B{2})
%driftwbw_1(driftwbw_1==0.05)=4;

% driftwbw_2 =size(driftwbw_1)
% driftwbw_2=log(driftwbw_1);



hp=responsive_stat & driftwbw_1==10;
lp= responsive_stat & driftwbw_1<0.5;

layerAgeCDF(driftwbw_1(:,1),age,lyr,inh,SF_pref_stat,'SF_pref tuning width Stationary');
% layerAgeCDF(drift_theta_w(:,2),age,lyr,inh,tuned_OS_w_run ,'OS tuning width Running');

layerAgePlot(driftwbw_1,age,lyr,inh,SF_pref_stat,'SF_pref tuning width Stationary');
% layerAgePlot(drift_theta_w(:,2),age,lyr,inh,tuned_OS_w_run,'OS tuning width Running');
p_wpref_bw=anovan(driftwbw_1(peak(:,1)>=2),g, 'interaction');


SF_pref_bw=SF_pref_stat & driftwbw_1<8;
layer_age_line(driftwbw_1,age,lyr,inh,SF_pref_stat,'SF_pref BW');
layerAgePlot(driftwbw_1,age,lyr,inh,SF_pref_stat,'SF_pref tuning width Stationary');
p_wpref_bw=anovan(driftwbw_1(lp),g, 'interaction');

layerAgePlot_SF_BW(driftwbw_1,age,lyr,inh,responsive_stat,{'pref Orient' 'Prct total'},'SF_pref tuning width Stationary');


%layerAgeScatterMedian(driftwbw_1(:,1),age,lyr,inh,bw_stat,'SF_pref tuning width Stationary');
% layerAgeScatterMedian(drift_theta_w(:,2),age,lyr,inh,tuned_OS_w_run,'OS tuning width Running');




%%%                         Receptive field STA data

%%%Nx data
clear f x f2 x1
good_STA =STA_exp_var>=0.6 & STA_ny<0.39;
responsive_stat = peak(:,1)>=2;  % firing rate (responsiveness) criteria for whether cells enter subsequent statistical analysis

clear STA_nx_EO1 STA_nx_adult
STA_nx_EO1=abs(STA_nx(age==1   & ~inh' & good_STA));
STA_nx_adult=abs(STA_nx(age==4  & ~inh' & good_STA));

n1=sum(~isnan(STA_nx_EO1))
n2=sum(~isnan(STA_nx_adult))
[h p]=vartest2(STA_nx_EO1,STA_nx_adult);

STA_ny_EO1=abs(STA_ny(age==1   & ~inh' & good_STA));
STA_ny_adult=abs(STA_ny(age==4  & ~inh' & good_STA));

[p h] = kstest2(STA_ny_EO1,STA_ny_adult)

clear p h g f g_s age_peak lyr_peak x1 y1
[p h]=lillietest(STA_nx_adult)
g_s=skewness(STA_nx_adult)
n=size(STA_nx_adult)

[p h]=lillietest(STA_nx_EO1)
g_s=skewness(STA_nx_EO1)
n=size(STA_nx_EO1)
%%calculate the binomial probability
 
%group1
clear total resp frac errdata prct_err pci phat
total=sum(age'==1  & lyr<=4  & ~inh & responsive_stat);
resp = sum(~isnan(STA_nx_EO1));         
frac= resp/total;
        
[phat,pci]= binofit(resp,total);        
% errdata =pci/sqrt(total);
% prct_err= errdata/resp;         
% prct_err_lin=prct_err*frac;
 %group 2
 clear total_2 resp_2 frac_2 errdata_2 prct_err_2 M2 V2
total_2=sum(age'==4  & lyr<=4  & ~inh & responsive_stat);
resp_2 = sum(~isnan(STA_nx_adult));         
frac_2= resp_2/total_2;
        
[M2,V2]= binofit(resp_2,total_2);        



%%ranksum
 %%%KS test
 B={STA_ny(age==1 & lyr'==4 & ~inh' & good_STA),STA_ny(age==4 & lyr'==4 & ~inh' & good_STA)};
[h p]= kstest2(B{1},B{2})

B={STA_nx(age==1 & lyr'==4 & ~inh' & good_STA),STA_nx(age==4 & lyr'==4 & ~inh' & good_STA)};
[h p]= kstest2(B{1},B{2})

n1=size(STA_nx(age==1 & lyr'<=6 & ~inh' & good_STA))
n2=size(STA_nx(age==4 & lyr'<=6 & ~inh' & good_STA))
%%%all layers
figure
[f,x]=hist(abs(STA_nx(age==1 & lyr'<=6 & ~inh' & good_STA )));
H1=bar(x,f/sum(f),'b');
hold on
[f1,x1]=hist(abs(STA_nx(age==4 & lyr'<=6 & ~inh'& good_STA )));
H2=bar(x,f1/sum(f1),'g');
ch=get(H2,'child');
set(ch,'facea',.5)
title 'Eo1_vs adult _nx'

figure
[f,x]=hist(abs(STA_ny(age==1 & lyr'<=6 & ~inh' & good_STA )));
H1=bar(x,f/sum(f),'b');
hold on
[f1,x1]=hist(abs(STA_ny(age==4 & lyr'<=6 & ~inh'& good_STA )));
H2=bar(x,f1/sum(f1),'g');
ch=get(H2,'child');
set(ch,'facea',.5)
title 'Eo1_vs adult _ny'
ranksum(abs(STA_ny(age==4 & lyr'<=6 & ~inh' & good_STA )),abs(STA_ny(age==2 & lyr'<=6 & ~inh'& good_STA )))


STA_nx=STA_nx';
STA_ny=abs(STA_ny)';
good_STA =STA_exp_var>=0.6 & STA_ny<0.39;

layerAgePlot(STA_nx,age,lyr,inh,good_STA','SF_pref tuning width Stationary nx');
layerAgePlot(abs(STA_ny),age,lyr,inh,good_STA','SF_pref tuning width Stationary ny');

jitterValuesX = 2*(rand(size(STA_nx))-0.5)*0.005;
STA_nx_j= STA_nx + jitterValuesX
STA_nx_j=STA_nx_j'

layerAgeNxNy(abs(STA_nx),abs(STA_ny),age,lyr,inh,good_STA,{'nx','ny'},'nx vs ny');
layerAgePlot_ratio_jlh(abs(STA_nx_j),abs(STA_ny),age,lyr,inh,good_STA',{'nx','ny'},'nx vs ny');


% layer_age_line(STA_nx(:,1),age,lyr,inh,good_STA','nxStationary ');
% layer_age_line(STA_ny(:,1),age,lyr,inh,good_STA','nxStationary ');


x1=STA_x_ratio(age==1 & lyr'==4 & good_STA & STA_ny'<0.39);
x2=STA_x_ratio(age==4 & lyr'==4 & good_STA & STA_ny'<0.39);
[p h]= ranksum(x1,x2);

x1=STA_x_ratio(age==1 & lyr'<=6 & ~inh' & good_STA & STA_ny'<0.39);
x2=STA_x_ratio(age==4 & lyr'<=6 & ~inh' & good_STA & STA_ny'<0.39);
[p h]= ranksum(x1,x2);

clear f f1 x x1
figure
[f,x]=hist(STA_ny(age==1   & ~inh' & good_STA ),0:0.035:0.4)
[f1,x1]=hist(STA_ny(age==4 & ~inh' & good_STA ),0:0.035:0.4)
bw=[(f/sum(f));(f1/sum(f1))];
figure
bar(x,bw',2)



% title 'adult'


%%% receptive field size from STAs

%%%all layers

good_STA =STA_exp_var>=0.60;
clear f x f1 x1

% A=(0.5).*STA_sigx;
% B=(0.5).*STA_sigy;
% area=pi.*A.*B;

A= 0.7031.*STA_sigx; %%0.7031 is the calculation for deg per pix 
B= 0.7031.*STA_sigy;

area=pi.*A.*B;
figure
hist(area(area<=100))



% med_EO1=nanmedian(area(age==1 & lyr'<=6  & good_STA & area<=100));
% s_EO1= semedian(area(age==1 & lyr'<=6  & good_STA & area<=100 ));
% 
% 
% med_adult=nanmedian(area(age==2& lyr'<=6  & good_STA & area<=100 ));
% s_adult= semedian(area(age==2 & lyr'<=6  & good_STA & area<=100));
% 
% figure
% barweb([med_EO1;med_adult],[s_EO1;s_adult]);
% title 'median area across all layers'
% [p,h]=ranksum(area(age==1 & lyr'==4  & good_STA & area<=100),area(age==2 & lyr'==4  & good_STA & area<=100));

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




x=nanmedian(A(age==1  &lyr'==4  & good_STA));
sx=semedian(A(age==1& lyr'==4 & good_STA))


x1=nanmedian(A(age==2 & lyr'==4 & good_STA));
sx1=semedian(A(age==2& lyr'==4 & good_STA))

figure
barweb([x;x1],[sx;sx1])

[p,h]=ranksum((A(age==1 & lyr'==4  & good_STA & area<=100)),(A(age==2 & lyr'==4  & good_STA & area<=100)))




x_B=nanmedian(B(age==1 & lyr'==4  & good_STA));
sx_B=semedian(B(age==1&lyr'==4  & good_STA))

x1_B=nanmedian(B(age==2 & lyr'==4  & good_STA));
sx1_B=semedian(B(age==2& lyr'==4 & good_STA))

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
layerAgePlot(A,age,lyr,inh,good_STA,'RF area');
layerAgePlot(B,age,lyr,inh,good_STA,'RF area');



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
tunedOSI_F1F0_stat =peak(:,1)>=1 & OSI(:,1)>=0.5 & driftF1F0(:,1) <2.2 ;
F1F0_stat = peak(:,1)>=2 & driftF1F0(:,1) <2.2;

% tunedOSI_F1F0_run = tunedOSI_run & driftF1F0(:,2) <=2.1 ;
 %F1F0_run = responsive_run & driftF1F0(:,2) <=2.1;

 F1= driftF1F0(:,1)>1;
layerAgeCDF(driftF1F0(:,1),age,lyr,inh,F1F0_stat,'F1F0 Stationary');% F1F ratio of all units regardless of selectivity
%layerAgeCDF(driftF1F0(:,2),age,lyr,inh,F1F0_run,'F1F0 Running');
 
layerAgePlot(driftF1F0(:,1),age,lyr,inh, F1F0_stat,'F1F0 Stationary');
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



