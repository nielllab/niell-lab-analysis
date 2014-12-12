%function compile developmental data
clear all
close all
dbstop if error
% dbclear all % this will exit out of db mode.
%[fname pname] =uiputfile('*.ps','pdf output'); psfilename=fullfile(pname,fname);  %%% get ps filename
psfilename = 'c:/test.ps';   %%% default location
if exist(psfilename,'file')==2;delete(psfilename);end %%% check for previous file

apath = 'E:\Jennifer_analysis\'; %apath = 'D:\Jen_ephys_data\developmental_periods\';
N =0; cells=0;  all_img_STA={};PINPed=0;

for dataset = 1:3  %%% control ("wt") NR5A1-cre/CHR2 animals vs. NR2A deleted NR5A1-cre
    
    if dataset ==1
        
        
        afiles = { 'NR2A\WT\8_11_14_wt_control\analysis_8_11_14_wt.mat',...
            'NR2A\WT\9_22_14_NR2A_wt\rec1\analysis_9_22_14_rec1A_wt.mat',...
            'NR2A\WT\10_9_14\analysis_10_9_14_wt_rec1.mat',...
            'NR2B\9_9_14_ctl\analysis_9_9_14_NR2B_wt_ctl.mat'};
        
    elseif dataset==2
        
        afiles = {'NR2A\KO\8_12_14_NR2A_KO\Rec2\analysis_8_12_14_rec2.mat',...
            'NR2A\KO\9_19_14\rec2\analysis_9_19_14_NR2A_rec2.mat',...
            'NR2A\KO\10_8_14\full clustering_rec1\analysis_10_8_14_NR2A_KO_rec1.mat'};
        
    elseif dataset==3
        afiles = {'NR2B\KO\9_6_14_NR2B_homo_mu\analysis_9_6_14_NR2B_homo_mu.mat',...
            'NR2B\KO\12_8_14_homo\analysis_12_8_14_NR2BKO.mat'};
        
    end
    
    
    
    for i = 1:length(afiles)
        
        clear params
        clear wn wn_movement
        clear LFP_movement
        clear bars
        clear wave_all
        clear rf_width
        clear locomotion
        
        load([apath afiles{i}]);
        
        clusterfilename
       
        afiles{i}
    
        if exist(clusterfilename,'file')
            clusterFile = clusterfilename;
        elseif exist(clusterfilename((length(pname)+1):end),'file')
            clusterFile = clusterfilename((length(pname)+1):end);
        elseif exist([clusterfilename((length(pname)+1):end) '.mat'],'file')
            clusterFile = [clusterfilename((length(pname)+1):end) '.mat'];
        else
            [fname pname] = uigetfile('*.mat','cluster file');
            clusterFile = fullfile(pname,fname);
            clusterfilename = clusterFile;
            if fname~=0
                save([apath afiles{i}],'clusterfilename','-append');
            end
        end
        clusterFile
        try
            load(clusterFile,'wave_all');
        catch
            display('no cluster file')
        end
        
        
        
        n_units = length(L_ratio);
        cellrange = N+1:N+n_units;
        N=N+n_units;
       
        
        number(i) = n_units;
        
        alldata( cellrange,1:2) = cells;
        alldata( cellrange,3) = L_ratio;
        
        pinp(cellrange,:)=PINPed';
        
        %%% waveform
        alldata( cellrange,4) = trough_width;
        alldata( cellrange,5) = trough2peak;
        alldata( cellrange,6) = -trough_depth./peak_height;
        alldata( cellrange,7:25)= wv';
        
     
        
        GT(cellrange)=4-dataset;
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
        driftpeak(cellrange,:)=field2array(drift,'maxFR');
        driftstd(cellrange,:)=field2array(drift,'peakstd');
        driftSigNoise(cellrange,:)=field2array(drift,'signoise');
         
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
        cvDSI(cellrange,:) = field2array(drift,'cv_dsi'); %%also circular variance measure change to "cv_dsi and cv_osi" in new compile programs
        cvOSI(cellrange,:)=field2array(drift,'cv_osi');
        %driftOri(cellrange,:) = field2array(drift,'orientfreq_all');
        
        if exist('bars');

        bar_spont(cellrange,:)=field2array(bars,'spont');
       
        else
        bar_spont(cellrange,:)= NaN;
        
        end
        
     clear meanwaves snr stdwaves
%         
        for c=1:length(cells);
            %%% get SNR
            tet =ceil(cells(c,1)/4);
            
            if exist('wave_all','var')
                wvall = wave_all{tet};
                wvclust = wvall(find(idx_all{(tet-1)*4+1}==cells(c,2)),:,:);
                
                amps =squeeze(min(wvclust(:,5:10,:),[],2));
                mn = abs(nanmean(amps));
                stdev = nanstd(amps);
                [y ind] = max(mn);
                snr(c) = mn(ind)/stdev(ind);
                
                meanwaves(c,:,:) = squeeze(nanmean(wvclust,1));
                stdwaves(c,:,:) = squeeze(nanstd(wvclust,[],1));
            else
                meanwaves=NaN;
                snr=NaN
                stdwaves=NaN
            end
            
%             
        end
        
        SNRall(cellrange)=snr;
        meanWavesAll(cellrange,:,:) = meanwaves;
        stdWavesAll(cellrange,:,:) = stdwaves;
        
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
       
        %% get peristimulus histograms for PINPed units
        psth_pinp(cellrange,:)=psth;
        histbins=histbins
     
        %get firing rate at all measured orients and SF, put into an array:
        %12 rows(orientations) by 7 columns(SpatialFreqs) for each cell
        
        drift_Ori_Sf(cellrange,:) = arrayfun(@(x)(getfield(x,'orientfreq_all')),drift,'UniformOutput',false);
       % drift_all(cellrange,:)=drift';
        
        if exist('locomotion');
           Vel(i,dataset)=arrayfun(@(x)(getfield(x,'mouseV')),locomotion,'UniformOutput',false);  
        end
       
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

figure
plot(alldata(GT==1,5),alldata(GT==1,6),'bo'); hold on
plot(alldata(GT==2,5),alldata(GT==2,6),'mo'); hold on
plot(alldata(GT==3,5),alldata(GT==3,6),'ro'); hold on
%plot(alldata(GT==4,5),alldata(GT==4,6),'go'); hold on

legend('EO','adult');

% [coeff score latent] = princomp(wvform);
% figure
% plot(latent);
% figure
% plot(score(:,1),score(:,2),'o');



inh = alldata(:,6)>0 & alldata(:,6)<4  & alldata(:,5)<8.5;  %%% directly based on wvform; k-means includes another inh group?
midnarrow = alldata(:,6)>0 & alldata(:,6)<4 & alldata(:,5)<12 &  alldata(:,5)>8.5;  %%% could analyze these specifically at some point
exc= alldata(:,5)>10 & alldata(:,6)>0 & alldata(:,6)<10;

figure
plot(alldata(find(inh),5),alldata(find(inh),6),'ro');
hold on
plot(alldata(find(exc),5),alldata(find(exc),6),'ko');
hold on
plot(alldata(find(pinp),5),alldata(find(pinp),6),'go');

figure
plot(wvform(exc,:)','color','k');hold on
plot(wvform(inh,:)','color','r');hold on
plot(wvform(pinp,:)','color','b');

for i = 1:size(driftA1,1)
    for j=1:size(driftA1,2)
        driftA1(i,j);
        driftA2(i,j);
        driftB(i,j);
        drift_theta_w(i,j);
        [OSI(i,j) DSI(i,j) width(i,j) peak(i,j)] = calculate_tuning(driftA1(i,j),driftA2(i,j),driftB(i,j),drift_theta_w(i,j));
        
        
    end
end


%variables created to sift through data conditionally

%%firing rate
responsive_run = peak(:,2)>=1; 
responsive_stat = peak(:,1)>=1;  % firing rate (responsiveness) criteria for whether cells enter subsequent statistical analysis

good_units=SNRall(:,1)>=2;

 pinp1= double(pinp);
 pinp1(pinp1==1)=2;
 pinp1(pinp1==0)=1;
 
figure
plot(histbins,psth_pinp(pinp,:));hold on
 
evoked_stat= peak(:,1)>=2 & good_units;  %OSI(:,1)>0.5 
evoked_run= responsive_run 
responsive_both= peak(:,1)>=1  & peak(:,2)>=1;

layerGTpinpPlot(peak(:,1),pinp1,GT,lyr,inh,good_units,'evoked FR');
layer_GT_line_pinp(peak(:,1),GT,lyr,inh,exc,pinp,good_units,'evoked FR');

j=pinp1==2 & GT'==2 & lyr==4 & exc & good_units

%%%ratio of running peak to stat peak
layerGenoPlot_pinp_ratio(peak(:,1),peak(:,2),pinp1,GT,lyr,inh,exc,evoked_stat,'gain-evokedFR');
layer_GT_pinp_bar_ratio(peak(:,1),peak(:,2),GT,lyr,inh,exc,pinp,evoked_stat,'gain-evokedFR' );


layerGTpinpPlot(driftspont(:,1),pinp1,GT,lyr,inh,exc,evoked_stat,'spont FR');
layerGTpinpPlot(driftspont(:,2),pinp1,GT,lyr,inh,exc,evoked_stat,'spont FR');

layer_GT_line_pinp(driftspont(:,1),GT,lyr,inh,exc,pinp,evoked_stat,'spont FR');
layer_GT_pinp_bar_ratio(driftspont(:,1),driftspont(:,2),GT,lyr,inh,exc,pinp,evoked_stat,'gain-spontFR' );



% layerAgePlot(rfw(:,1),GT,lyr,inh,responsive_stat,'RF size');
% layer_age_line(rfw(:,1),GT,lyr,inh,responsive_stat,'RF size');
 
x1=responsive_stat & ~inh & GT'==1 & driftspont1(:,1)>=0.045 & lyr<=3;
y1=responsive_stat & ~ inh & GT'==4 & driftspont1(:,1)>=0.045 &lyr<=3;
[p h]=ranksum(driftspont1(x1),driftspont1(y1));

%%preferred oprientation 
tunedOSI_stat = evoked_stat & OSI (:,1)>=0.5; 
drift_theta_1=(drift_theta*180)/pi;
drift_theta_1(drift_theta_1>330)=0;
layerGT_pinp_Plot_pref_Orient(drift_theta_1(:,1),pinp1,GT,lyr,inh,exc,tunedOSI_stat ,'Prefered Orientation Stationary');

%%%OS
%layerAgeCDF(OSI(:,1),GT,lyr,inh,responsive_stat,'OSI all Stationary');
layerGTpinpPlot(OSI(:,1),pinp1,GT,lyr,inh,exc,evoked_stat,'OSI Stationary');
layer_GT_line_pinp(OSI(:,1),GT,lyr,inh,exc,pinp,evoked_stat,'spont FR');

%%
% %%%DS
layerGTpinpPlot(cvDSI(:,1),pinp1,GT,lyr,inh,exc,evoked_stat,'DSI Stationary');

layer_GT_line(cvDSI(:,1),GT,lyr,inh,DSI_stat_top50,'DSI top 50%');
layer_GT_line(cvDSI(:,1),GT,lyr,inh,DSI_stat,'DSI ');

%%%                              simple cells F1F0
tunedOSI_F1F0_stat =evoked_stat & OSI(:,1)>=0.5 ;

layerGTpinpPlot(driftF1F0(:,1),pinp1,GT,lyr,inh, exc,tunedOSI_F1F0_stat,'F1F0 Stationary');
%layerAgePlot(driftF1F0(:,2),age,lyr,inh, F1F0_run,'F1F0 Running');

layerGT_pinp_Plot_frac_simple(driftF1F0(:,1),pinp1,GT,lyr,inh,exc, tunedOSI_F1F0_stat,'F1F0 Stationary');
%layerAgePlot_frac_simple(driftF1F0(:,2),age,lyr,inh, F1F0_run,'F1F0 Running');


%%%                         Receptive field STA data

%%%Nx data
clear f x f2 x1
good_STA =STA_exp_var>=0.6 & STA_ny<0.39;
responsive_stat = peak(:,1)>=2;  % firing rate (responsiveness) criteria for whether cells enter subsequent statistical analysis







