%function compile DOI/drug treatment data
clear all
close all
dbstop if error
% dbclear all % this will exit out of db mode.
%[fname pname] =uiputfile('*.ps','pdf output'); psfilename=fullfile(pname,fname);  %%% get ps filename
psfilename = 'c:/test.ps';   %%% default location
if exist(psfilename,'file')==2;delete(psfilename);end %%% check for previous file

apath = 'D:\'; %apath = 'D:\Jen_ephys_data\developmental_periods\';
N =0; cells=0;  all_img_STA={};PINPed=0; STA_peak={};stopCRF={}; moveCRF={};
sessionNum=0;
clear wn
psfilename = 'D:\Angie_analysis\analysisPS.ps';
if exist(psfilename,'file')==2;delete(psfilename);end

for dataset = 1:1  %%% pre and post DOI
    
    if dataset==1
        %%DOI
%         
         afiles = {'Angie_analysis\DOI_experiments\08_31_15\analysis_083115.mat'} %,...
% %             'Angie_analysis\DOI_experiments\10_28_15\analysis_102815.mat'}
% %         %    'Angie_analysis\DOI_experiments\11_04_15\analysis_110415.mat'

    
% %         
%         afiles = {'Angie_analysis\DOI_experiments\7_29_15\analysis_072915.mat',...
%             'Angie_analysis\DOI_experiments\08_02_15\analysis_080215.mat',...
%             'Angie_analysis\DOI_experiments\08_03_15\analysis_080315.mat',...
%             'Angie_analysis\DOI_experiments\08_05_15\analysis_080515.mat',...
%             'Angie_analysis\DOI_experiments\08_07_15\analysis_080715.mat',...
%             'Angie_analysis\DOI_experiments\09_02_15\analysis_090215.mat',...
%             'Angie_analysis\DOI_experiments\09_04_15\analysis_090415.mat',...
%             'Angie_analysis\DOI_experiments\09_07_15\analysis_090715.mat',...
%             'Angie_analysis\DOI_experiments\09_08_15\analysis_090815.mat',...
%             'Angie_analysis\DOI_experiments\09_23_15\analysis_092315.mat',...
%             'Angie_analysis\DOI_experiments\09_29_15\analysis_092915.mat',...
%             'Angie_analysis\DOI_experiments\09_30_15\analysis_093015.mat',...
%             'Angie_analysis\DOI_experiments\10_19_15\analysis_101915.mat',...
%             'Angie_analysis\DOI_experiments\10_28_15\analysis_102815.mat',...
%             'Angie_analysis\DOI_experiments\11_04_15\analysis_110415.mat',...
%             'Angie_analysis\DOI_experiments\11_19_15\analysis_111915.mat'}
%      
        
        
        % %                 %lisuride
        
        %               afiles = {'Angie_analysis\DOI_experiments\08_13_15\08_13_15a\analysis_081315a.mat',...
        % %                     'Angie_analysis\DOI_experiments\08_25_15\analysis_082515.mat',...
        % %                     'Angie_analysis\DOI_experiments\08_26_15\analysis_082615.mat',...
        % %                     'Angie_analysis\DOI_experiments\08_31_15\analysis_083115.mat'}
        %
        % saline
%         
%         afiles = {'Angie_analysis\DOI_experiments\08_12_15\analysis_081215.mat',...
%             'Angie_analysis\DOI_experiments\08_21_15\analysis_082115.mat',...
%             'Angie_analysis\DOI_experiments\09_14_15\analysis_091415.mat',...
%             'Angie_analysis\DOI_experiments\09_18_15\analysis_091815.mat',...
%             'Angie_analysis\DOI_experiments\09_21_15\analysis_092115.mat'}
%    
    elseif dataset==2
        %         %%lisuride
        %afiles = {'Angie_analysis\DOI_experiments\08_13_15\08_13_15a\analysis_081315a.mat',...
        %             'Angie_analysis\DOI_experiments\08_25_15\analysis_082515.mat',...
        %             'Angie_analysis\DOI_experiments\08_26_15\analysis_082615.mat',...
        %             'Angie_analysis\DOI_experiments\08_31_15\analysis_083115.mat'}
        %
    elseif dataset==3
        
        %%EtOH\DMSO\Saline
       % afiles = {'Angie_analysis\DOI_experiments\7_16_15\analysis_etoh_control.mat'};
        
        
        
    end
    
    
    
    for i = 1:length(afiles)
        
        sessionNum=sessionNum+1;
        clear params
        clear wn
        clear wn_post
        
        clear wn_movement
        clear LFP_movement
        clear bars
        clear wave_all
        clear rf_width
        clear locomotion
        clear drift
        clear layer
        clear psth
        clear psth_power2
        
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
        
        session(cellrange)=sessionNum;
        number(i) = n_units;
        
        alldata( cellrange,1:2) = cells;
        alldata( cellrange,3) = L_ratio;
        
        %pinp(cellrange,:)=PINPed';
        
        %%% waveform
        alldata( cellrange,4) = trough_width;
        alldata( cellrange,5) = trough2peak;
        alldata( cellrange,6) = -trough_depth./peak_height;
     %   alldata( cellrange,7:25)= wv';
        alldata(cellrange,26)=layer;
        
        
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
        
        %         if exist('rf_width');
        %
        %         rfw(cellrange,:) = rf_width*30;
        %
        %         else
        %
        %         rfw(cellrange,:) = NaN;
        %         end
        
        
        if  exist('drift', 'var');
            
            driftA1(cellrange,:)= field2array(drift,'A1');
            driftA2(cellrange,:)=field2array(drift,'A2');
            driftB(cellrange,:)= field2array(drift,'B');
            %driftpeak(cellrange,:)=field2array(drift,'maxFR');
            %driftstd(cellrange,:)=field2array(drift,'peakstd');
            %driftSigNoise(cellrange,:)=field2array(drift,'signoise');
            
            drift_theta_w(cellrange,:)=field2array(drift,'thetawidth');
            drift_theta(cellrange,:)=field2array(drift,'theta');
            
            driftspont(cellrange,:) = field2array(drift,'spont');
            
            driftwpref(cellrange,:) = field2array(drift,'wpref');
            driftwbw(cellrange,:) = field2array(drift,'bw') ;
            
            driftF1F0(cellrange,:) = field2array(drift,'F1')./field2array(drift,'F0');
            driftF0(cellrange,:) = field2array(drift,'F0');
            %       driftorientfreq_all(cellrange,:)=field2array(drift, 'orientfreq_all');
            
%             cvDSI(cellrange,:) = field2array(drift,'cv_dsi'); %%also circular variance measure change to "cv_dsi and cv_osi" in new compile programs
%             cvOSI(cellrange,:)=field2array(drift,'cv_osi');
            %driftOri(cellrange,:) = field2array(drift,'orientfreq_all');
        else
            display('no drift???')
            driftA1(cellrange,1:2)= NaN;
            driftA2(cellrange,1:2)=NaN;
            driftB(cellrange,1:2)= NaN;
%             driftpeak(cellrange,1:2)=NaN;
%             driftstd(cellrange,1:2)=NaN;
%             driftSigNoise(cellrange,1:2)=NaN;
            
            drift_theta_w(cellrange,1:2)=NaN;
            drift_theta(cellrange,1:2)=NaN;
            
            driftspont(cellrange,1:2) = NaN;
            
            driftwpref(cellrange,1:2) = NaN;
            driftwbw(cellrange,1:2) = NaN ;
            
            driftF1F0(cellrange,1:2) = NaN;
            driftF0(cellrange,1:2) = NaN;
%            driftorientfreq_all(cellrange,:)=field2array(drift, 'orientfreq_all');
%             cvDSI(cellrange,1:2) = NaN; %%also circular variance measure change to "cv_dsi and cv_osi" in new compile programs
%             cvOSI(cellrange,1:2)=NaN;
        end
% %         
        if exist('bars','var'); %#ok<EXIST>
            
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
            % %
            % % %
        end
        %
        SNRall(cellrange)=snr;
        meanWavesAll(cellrange,:,:) = meanwaves;
        stdWavesAll(cellrange,:,:) = stdwaves;
        
        if exist('params','var');
            
            
            all_img_STA(cellrange)= all_img;
            all_fit_STA(cellrange)=all_fit;
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
        
        %         if exist ('wn_movement','var')
        %
        %             stopCRF(cellrange)=field2array_false(wn_movement,'stopCRF');
        %             moveCRF(cellrange)=field2array_false(wn_movement,'moveCRF');
        %         else
        %             stopCRF(cellrange)=NaN;
        %             moveCRF(cellrange)=NaN;
        %         end
        
        clear m ind x y t_lag STA STA1 STA_post STA1_post c crf crf1 crf_post crf1_post
        
        if exist ('wn','var')
            for w = 1:length(wn)
                STA = wn(w).sta;
                
                %%%Dtermine time point with maximial response
                [m ind] = max(abs(STA(:)-127));
                [x y t_lag] = ind2sub(size(STA),ind);
                
                STA1{w} = STA(:,:,t_lag)-128;
                % figure
                % imagesc(STA1{1,w}',[-64 64]); axis equal
            end
            STA_peak(cellrange)=STA1
            
            for c = 1:length(wn)
                if isempty(wn(c).crf)
                    wn(c).crf = zeros(20,1);
                    wn(c).spont = 0;
                end
                crf = wn(c).crf;
                crf1{c} = crf;
            end
            
            CRF(cellrange)=crf1;
            wn_spont(cellrange,:)= field2array_false(wn,'spont');
            
        else
            STA_peak(cellrange)=NaN
            CRF(cellrange)=NaN
            wn_spont(cellrange)=NaN
        end
        %close all
        
        clear w c
        
        if exist ('wn_post','var')
            for w = 1:length(wn_post)
                STA_post = wn_post(w).sta;
                
                %%%Dtermine time point with maximial response
                [m ind] = max(abs(STA_post(:)-127));
                [x y t_lag] = ind2sub(size(STA_post),ind);
                
                STA1_post{w} = STA_post(:,:,t_lag)-128;
            end
            
            for c = 1:length(wn_post)
                if isempty(wn_post(c).crf)
                    wn_post(c).crf = zeros(20,1);
                    wn_post(c).spont = 0;
                end
                crf_post = wn_post(c).crf;
                crf1_post{c} = crf_post;
            end
            
            CRF_post(cellrange)=crf1_post;
            STA_peak_post(cellrange)=STA1_post
            wn_spont_post(cellrange,:)= field2array_false(wn_post,'spont');
        else
            STA_peak_post(cellrange)=NaN
            wn_spont_post(cellrange)=NaN
        end
        %size(wvform)
        wvform(cellrange,:) = wv';
        
        
        
    end %%% loop over adult vs EO
end

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


inh = alldata(:,6)>0 & alldata(:,6)<1.6  & alldata(:,5)<7.5  %%% directly based on wvform; k-means includes another inh group?
midnarrow = alldata(:,6)>0 & alldata(:,6)<4 & alldata(:,5)<10 &  alldata(:,5)>7.5;  %%% could analyze these specifically at some point
exc= alldata(:,5)>10 & alldata(:,6)>0 & alldata(:,6)<4;

inh = inh | midnarrow;
lyr = alldata(:,26)

figure
subplot(2,3,1)
plot(alldata(find(inh),5),alldata(find(inh),6),'ro');
hold on
plot(alldata(find(exc),5),alldata(find(exc),6),'ko');
plot(alldata(find(midnarrow),5),alldata(find(midnarrow),6),'bo');
xlim([5 16])
ylim([1 3.5])
title('all units')

for selectlyr=3:6
    
    subplot(2,3,selectlyr-1)
    plot(alldata(find(inh & lyr==selectlyr),5),alldata(find(inh & lyr==selectlyr),6),'ro');
    hold on
    plot(alldata(find(exc &lyr==selectlyr),5),alldata(find(exc & lyr==selectlyr),6),'ko');
    plot(alldata(find(midnarrow&lyr==selectlyr),5),alldata(find(midnarrow & lyr==selectlyr),6),'bo');
    xlim([5 16])
    ylim([1 3.5])
    title(sprintf('layer %d',selectlyr))
end


% OSI = zeros(size(driftA1)); OSI(:) = NaN;
% DSI = OSI; width=OSI; peak=OSI;
for i = 1:size(driftA1,1)
    if ~isnan(driftA1(i,1))
        for j=1:size(driftA1,2)
        driftA1(i,j);
        driftA2(i,j);
        driftB(i,j);
        drift_theta_w(i,j);
        [OSI(i,j) DSI(i,j) width(i,j) peak(i,j)] = calculate_tuning(driftA1(i,j),driftA2(i,j),driftB(i,j),drift_theta_w(i,j));


        end
    end
 end

%%% define groups
DOI = GT'==3;
Lis=GT'==2;
EtOH = GT'==1;



%%%plotting CRF for units
%meancontrast = 0.5*(contrastdata(1:10)+contrastdata(20:-1:11));

use = find(DOI)
clear crf
figure
for i = 1:length(use)
    subplot(26,25,i); hold on
    set(gcf,'Position',[300 300 800 400])
    
    
    meandata = 0.5*(CRF{:,use(i)}(1:10)+CRF{:,use(i)}(20:-1:11));
    if inh(use(i))
        plot(meandata,'r:');
    else
        plot(meandata,'b:');
    end
    crf(i,:,1) = meandata;
    meandata2 = 0.5*(CRF_post{:,use(i)}(1:10)+CRF_post{:,use(i)}(20:-1:11));
    if inh(use(i))
        plot(meandata2,'r');
    else
        plot(meandata2,'b');
    end
    ylim([0 max(5,1.2*max(max(meandata,meandata2)))])
    crf(i,:,2) = meandata2;
end
crf(8,:,:) = crf(8,:,:) +crf(9,:,:); crf(9,:,:)=NaN;
crf(10,:,:)= crf(10,:,:) + crf(11,:,:); crf(11,:,:)=NaN;


%spontaneous activity by layer
clear spont sponterr
figure
subplot(2,3,1)
spont(2,:) = squeeze(nanmean(crf(inh,1,:),1))
sponterr(2,:) = squeeze(nanstd(crf(inh,1,:),1))/sqrt(sum(inh));
spont(1,:) = squeeze(nanmean(crf(~inh,1,:),1))
sponterr(1,:) = squeeze(nanstd(crf(~inh,1,:),1))/sqrt(sum(~inh));
barweb(spont,sponterr)
set(gca,'Xticklabel',{'exc','inh'}); legend({'pre','post'})
ylabel('spont sp/sec')
title('all units')

for selectlyr=3:6
subplot(2,3,selectlyr)
spont(2,:) = squeeze(nanmean(crf(inh & lyr==selectlyr,1,:),1))
sponterr(2,:) = squeeze(nanstd(crf(inh & lyr==selectlyr,1,:),1))/sqrt(sum(inh& lyr==selectlyr));
spont(1,:) = squeeze(nanmean(crf(~inh & lyr==selectlyr,1,:),1))
sponterr(1,:) = squeeze(nanstd(crf(~inh & lyr ==selectlyr,1,:),1))/sqrt(sum(~inh& lyr==selectlyr));
barweb(spont,sponterr)
set(gca,'Xticklabel',{'exc','inh'})
ylabel('spont sp/sec')
title(sprintf('lyr %d',selectlyr))
end


%evoked activity by layer
clear resp resperr
figure
subplot(2,3,1)
resp(2,:) = squeeze(nanmean(crf(inh,10,:),1))
resperr(2,:) = squeeze(nanstd(crf(inh,10,:),1))/sqrt(sum(inh));
resp(1,:) = squeeze(nanmean(crf(~inh,10,:),1))
resperr(1,:) = squeeze(nanstd(crf(~inh,10,:),1))/sqrt(sum(~inh));
barweb((resp-spont),resperr)
hold on
set(gca,'Xticklabel',{'exc','inh'}); legend({'pre','post'})
ylabel('evoked sp/sec')
title('all layers')

for selectlyr=3:6
resp(2,:) = squeeze(nanmean(crf(inh & lyr==selectlyr,10,:),1))
resperr(2,:) = squeeze(nanstd(crf(inh & lyr==selectlyr,10,:),1))/sqrt(sum(inh& lyr ==selectlyr));
resp(1,:) = squeeze(nanmean(crf(~inh & lyr==selectlyr,10,:),1))
resperr(1,:) = squeeze(nanstd(crf(~inh & lyr==selectlyr,10,:),1))/sqrt(sum(~inh& lyr ==selectlyr));
subplot(2,3,selectlyr)
barweb((resp-spont),resperr)
set(gca,'Xticklabel',{'exc','inh'})
ylabel('evoked sp/sec')
title(sprintf('lyr %d',selectlyr))
end


%%% evoked scatter plots
figure
subplot(2,3,1)
hold on
%plot(crf(exc & pinped,10,1)-crf(exc& pinped,1,1),crf(exc& pinped,10,2)-crf(exc& pinped,1,2),'g.');
plot(crf(exc,10,1)-crf(exc,1,1),crf(exc,10,2)-crf(exc,1,2),'b.');
plot([0 10],[0 10]); axis equal
title('all exc; evoked')

subplot(2,3,2)
plot(crf(inh,10,1)-crf(inh,1,1),crf(inh,10,2)-crf(inh,1,2),'r.');
hold on
% plot(crf(inh & pinped,10,1)-crf(inh& pinped,1,1),crf(inh& pinped,10,2)-crf(inh& pinped,1,2),'g.');
plot([0 10],[0 10]); axis equal
title('all inh')

for selectlyr=3:6
    subplot(2,3,selectlyr)
    plot(crf(exc&lyr==selectlyr,10,1) - crf(exc&lyr==selectlyr,1,1), crf(exc&lyr==selectlyr,10,2) - crf(exc&lyr==selectlyr,1,2),'b.')
    hold on
    plot([0 10],[0 10])
    axis equal
    title(sprintf('lyr %d',selectlyr))
end


%%% spont scatter plots
figure
subplot(2,3,1)
hold on
plot(crf(exc,1,1),crf(exc,1,2),'b.');
plot([0 10],[0 10]); axis equal
title('all exc; spont')

subplot(2,3,2)
plot(crf(inh,1,1),crf(inh,1,2),'r.');
hold on
plot([0 10],[0 10]); axis equal
title('all inh')

for selectlyr=3:6
    subplot(2,3,selectlyr)
    plot( crf(exc&lyr==selectlyr,1,1), crf(exc&lyr==selectlyr,1,2),'b.')
    hold on
    plot([0 5],[0 5])
    axis equal
    title(sprintf('lyr %d',selectlyr))
end

%%% drift spont scatter plots
figure
subplot(2,3,1)
hold on
plot(driftspont(exc,1),driftspont(exc,2),'b.');
plot([0 10],[0 10]); axis equal
title('all exc; drift spont')

subplot(2,3,2)
plot(driftspont(inh,1),driftspont(inh,2),'r.');
hold on
plot([0 10],[0 10]); axis equal
title('all inh')

for selectlyr=3:6
    subplot(2,3,selectlyr)
    plot( driftspont(exc&lyr==selectlyr,1), driftspont(exc&lyr==selectlyr,2),'b.')
    hold on
    plot([0 5],[0 5])
    axis equal
    title(sprintf('lyr %d',selectlyr))
end


%%% drift peak scatter plots
figure
subplot(2,3,1)
hold on
plot(peak(exc,1),peak(exc,2),'b.');
plot([0 10],[0 10]); axis equal
title('all exc; drift peak')

subplot(2,3,2)
plot(peak(inh,1),peak(inh,2),'r.');
hold on
plot([0 10],[0 10]); axis equal
title('all inh')

for selectlyr=3:6
    subplot(2,3,selectlyr)
    plot( peak(exc&lyr==selectlyr,1), peak(exc&lyr==selectlyr,2),'b.')
    hold on
    plot([0 5],[0 5])
    axis equal
    title(sprintf('lyr %d',selectlyr))
end


%%%Jen's version of plots below need to make compatible with Cris above

for i = 1:ceil(length((use))/8)
    
    figure
    set(gcf,'Position',[10 10 500 800]);
    
    for j= 1:min(8,length(use)-(i-1)*8)
        
        if inh(use((i-1)*8+j))
            col = 'r';
        else col = 'b';
        end
        
        subplot(8,3,3*(j-1)+1);
        plot(wvform(use((i-1)*8+j),:),col);axis off
        
        subplot(8,3,3*(j-1)+2);
        
        if ~isempty(STA_peak{1,use((i-1)*8+j)})
            colormap jet
            imagesc(STA_peak{1,use((i-1)*8+j)}',[-64 64]);
        end
        title('Pre DOI');
        set(gca,'Ytick',[]); set(gca,'Xtick',[])
        
        subplot(8,3,3*(j-1)+3);
        meandata = 0.5*(CRF{:,use((i-1)*8+j)}(1:10)+CRF{:,use((i-1)*8+j)}(20:-1:11));
        if inh(use(i))
            plot(meandata,'r:');
        else
            plot(meandata,'b:');
        end
        
    end
    %
end



% use_lis = find( resp & Lis);
% for i=1:ceil(length(use_lis)/8);
%       set(gcf,'Position',[10 50 400 600]);
%       for j= 1:min(8,length(use_lis)-(i-1)*8)
%
%     subplot(8,8,3*(j-1)+1);
%    if inh(use_lis((i-1)*8+j))
%        col = 'r';
%    else col = 'b';
%    end
%     plot(wvform(use_lis((i-1)*8+j),:),col);axis off
%
%     subplot(8,3,3*(j-1)+2);
%     %plot(plotrange,psth_pinp(use((i-1)*8+j),plotrange));ylim([0 50]); xlim([min(plotrange) max(plotrange)]);
%
%     if ~isempty(STA_peak{1,use_lis((i-1)*8+j)})
%     colormap jet
%     imagesc(STA_peak{1,use_lis((i-1)*8+j)}',[-64 64]);
%     end
%     title('pre Lis');
%     %   title(sprintf('%0.0f',cells(use((i-1)*8+j))));
%   set(gca,'Ytick',[]); set(gca,'Xtick',[])
%
%   subplot(8,3,3*(j-1)+3);
%   meandata = 0.5*(CRF{:,use_lis((i-1)*8+j)}(1:10)+CRF{:,use_lis((i-1)*8+j)}(20:-1:11));
%     if inh(use_lis(i))
%         plot(meandata,'r:');
%     else
%         plot(meandata,'b:');
%     end
%        end
% end


%%plot waveforms and STAs post drug treatment



%use = find( resp & DOI);
for i=1:ceil(length(use)/8);
    figure
    set(gcf,'Position',[10 10 500 800]);
    for j= 1:min(8,length(use)-(i-1)*8)
        
        subplot(8,3,3*(j-1)+1);
        if inh(use((i-1)*8+j))
            col = 'r';
        else col = 'b';
        end
        plot(wvform(use((i-1)*8+j),:),col);axis off
        
        subplot(8,3,3*(j-1)+2);
        
        %plot(plotrange,psth_pinp(use((i-1)*8+j),plotrange));ylim([0 50]); xlim([min(plotrange) max(plotrange)]);
        
        if ~isempty(STA_peak_post{1,use((i-1)*8+j)})
            colormap jet
            imagesc(STA_peak_post{1,use((i-1)*8+j)}',[-64 64]);
        end
        title('Post DOI');
        set(gca,'Ytick',[]); set(gca,'Xtick',[])
        subplot(8,3,3*(j-1)+3);
        
        
        meandata2 = 0.5*(CRF_post{:,use((i-1)*8+j)}(1:10)+CRF_post{:,use((i-1)*8+j)}(20:-1:11));
        if inh(use(i))
            plot(meandata2,'r');
        else
            plot(meandata2,'b');
        end
        ylim([0 max(5,1.2*max(max(meandata,meandata2)))])
        crf(i,:,2) = meandata2;
        
    end
    
end

% %use_lis = find( resp & Lis);
% for i=1:ceil(length(use_lis)/8);
%        figure
%        set(gcf,'Position',[10 50 400 600]);
%        for j= 1:min(8,length(use_lis)-(i-1)*8)
%
%     subplot(8,8,3*(j-1)+1);
%     plot(wvform(use_lis((i-1)*8+j),:));axis off
%
%     subplot(8,8,3*(j-1)+2);
%
%     %plot(plotrange,psth_pinp(use((i-1)*8+j),plotrange));ylim([0 50]); xlim([min(plotrange) max(plotrange)]);
%
%     if ~isempty(STA_peak_post{1,use_lis((i-1)*8+j)})
%     colormap jet
%     imagesc(STA_peak_post{1,use_lis((i-1)*8+j)}',[-64 64]);
%     end
%     title('post Lis');
%     %   title(sprintf('%0.0f',cells(use((i-1)*8+j))));
%   set(gca,'Ytick',[]); set(gca,'Xtick',[])
%        end
% end


% figure
% plot(wn_spont1(~inh),wn_spont1_post(~inh),'bo')
% hold on
% plot(wn_spont1(inh),wn_spont1_post(inh),'ro')
% plot([0 8],[0 8]);
% axis square
% axis([0 8 0 8])

% PlotSpontData(wn_spont1,wn_spont1_post,DOI,1,inh,resp,'DOI')
% PlotSpontData(wn_Spont1,wn_spont_post1,Lis,lyr,inh,resp,'Lis')

%%pre
% use = find( resp & DOI);
% for i=1:ceil(length(use)/8);
%        figure
%        for j= 1:8
%
%     subplot(8,8,3*(j-1)+1);
%     plot(wvform(use((i-1)*8+j),:));axis off
%
%     subplot(8,8,3*(j-1)+2);
%
% if ~isempty(CRF{:,use((i-1)*8+j)})
% meandata = 0.5*(CRF{:,use((i-1)*8+j)}(1:10)+CRF{:,use((i-1)*8+j)}(20:-1:11));
% plot(meandata,'g');
% xlabel('contrast');
% ylabel('response');
% hold on;
% plot(CRF{:,use((i-1)*8+j)});
% end
%
%   title('pre DOI');
%
%        end
% end


% use = find(resp & DOI);
% for i=1:ceil(length(use)/8);
%        figure
%        for j= 1:8
%
%     subplot(4,4,3*(j-1)+1);
%     plot(wvform(use((i-1)*8+j),:));axis off
%
%     subplot(4,4,3*(j-1)+2);
%
% if ~isempty(CRF{:,use((i-1)*8+j)})
% meandata = 0.5*(CRF{:,use((i-1)*8+j)}(1:10)+CRF{:,use((i-1)*8+j)}(20:-1:11));
% plot(meandata,'g');
% xlabel('contrast');
% ylabel('response');
% hold on;
% plot(CRF{:,use((i-1)*8+j)});
% end
%
%   title('pre DOI');
%       figure
% for i = 1:length(use)
%     subplot(4,7,i); hold on
%
%     meandata = 0.5*(CRF{:,use(i)}(1:10)+CRF{:,use(i)}(20:-1:11));
%     if inh(use(i))
%         plot(meandata,'r:');
%     else
%         plot(meandata,'b:');
%     end
%     crf(i,:,1) = meandata;
%     meandata2 = 0.5*(CRF_post{:,use(i)}(1:10)+CRF_post{:,use(i)}(20:-1:11));
%     if inh(use(i))
%         plot(meandata2,'r');
%     else
%         plot(meandata2,'b');
%     end
%     ylim([0 max(5,1.2*max(max(meandata,meandata2)))])
%     crf(i,:,2) = meandata2;
% end
%        end
% end


%%post
% use = find(resp & DOI);
% for i=1:ceil(length(use)/8);
%        figure
%        for j= 1:8
%
%     subplot(4,4,2*(j-1)+1);
%     plot(wvform(use((i-1)*8+j),:));axis off
%
%     subplot(4,4,2*(j-1)+2);
%
% if ~isempty(CRF_post{:,use((i-1)*8+j)})
% meandata = 0.5*(CRF_post{:,use((i-1)*8+j)}(1:10)+CRF_post{:,use((i-1)*8+j)}(20:-1:11));
% plot(meancontrast,meandata,'g');
% xlabel('contrast');
% ylabel('response');
% hold on;
% plot(contrastdata,CRF_post{:,use((i-1)*8+j)});
% end
%  title('post DOI');
%
%        end
% end


%%%liseride CRFs

% %%%pre
% use = find(resp & Lis);
% for i=1:ceil(length(use)/8);
%        figure
%        for j= 1:8
%
%     subplot(4,4,2*(j-1)+1);
%     plot(wvform(use((i-1)*8+j),:));axis off
%
%     subplot(4,4,2*(j-1)+2);
%
% if ~isempty(CRF{:,use((i-1)*8+j)})
% meandata = 0.5*(CRF{:,use((i-1)*8+j)}(1:10)+CRF{:,use((i-1)*8+j)}(20:-1:11));
% plot(meancontrast,meandata,'g');
% xlabel('contrast');
% ylabel('response');
% hold on;
% plot(contrastdata,CRF{:,use((i-1)*8+j)});
% end
%
%   title('pre Lis');
%        end
% end

% %%%post
% use = find(resp & Lis);
% for i=1:ceil(length(use)/8);
%        figure
%        for j= 1:8
%
%     subplot(4,4,2*(j-1)+1);
%     plot(wvform(use((i-1)*8+j),:));axis off
%
%     subplot(4,4,2*(j-1)+2);
%
% if ~isempty(CRF_post{:,use((i-1)*8+j)})
% meandata = 0.5*(CRF_post{:,use((i-1)*8+j)}(1:10)+CRF_post{:,use((i-1)*8+j)}(20:-1:11));
% plot(meancontrast,meandata,'g');
% xlabel('contrast');
% ylabel('response');
% hold on;
% plot(contrastdata,CRF_post{:,use((i-1)*8+j)});
% end
%
%   title('post Lis');
%   %set(gca,'Ytick',[]); set(gca,'Xtick',[])
%        end
% end

%%% plot grating respons data
% layerTreatmentPlot(peak(:,1),GT,lyr,inh,resp,'peak')


% %%% waveform of cell of interest
% figure
% plot(wvform(exc,:)','b');
% hold on
% plot(wvform(inh,:)','r'); hold on
% plot(wvform(DOI,:)','g')


%generate STA for cells with the criteria specified by l

% l=find(~inh & DOI & resp);
%
% for j=1:length(STA_peak)
%
%     if ismember(j,l)
% figure
% if ~isempty(STA_peak{1,j})
%     colormap jet
%     imagesc(STA_peak{1,j}',[-64 64]);
% end
%     end
% end


% driftF1F0(629,1)
% OSI(629,1)
% DSI(629,1)
% driftwpref(629,1)

