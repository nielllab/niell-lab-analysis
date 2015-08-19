%function compile developmental data
clear all
close all
dbstop if error
% dbclear all % this will exit out of db mode.
%[fname pname] =uiputfile('*.ps','pdf output'); psfilename=fullfile(pname,fname);  %%% get ps filename
psfilename = 'c:/test.ps';   %%% default location
if exist(psfilename,'file')==2;delete(psfilename);end %%% check for previous file

apath = 'D:\Jen_analysis\'; %apath = 'D:\Jen_ephys_data\developmental_periods\';
N =0; cells=0;  all_img_STA={};PINPed=0; STA_peak={};stopCRF={}; moveCRF={};
sessionNum=0;

for dataset = 1:2  %%% control ("wt") NR5A1-cre/CHR2 animals vs. NR2A deleted NR5A1-cre
    
    if dataset ==1
        %%DOI
        
        afiles = { 'DOI_next_gen\7_2_15\analysis.mat'};

    elseif dataset==2
%         %%lisuride 
        afiles = {'DOI_next_gen\7_14_15\analysis.mat'};
%         
    elseif dataset==3
        
        %%EtOH\DMSO\Saline
        afiles = {''};
%          
           
        
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
        alldata( cellrange,7:25)= wv';
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

        cvDSI(cellrange,:) = field2array(drift,'cv_dsi'); %%also circular variance measure change to "cv_dsi and cv_osi" in new compile programs
        cvOSI(cellrange,:)=field2array(drift,'cv_osi');
        %driftOri(cellrange,:) = field2array(drift,'orientfreq_all');
        else
        driftA1(cellrange,:)= NaN;
        driftA2(cellrange,:)=NaN;
        driftB(cellrange,:)= NaN;
        driftpeak(cellrange,:)=NaN;
        driftstd(cellrange,:)=NaN;
        driftSigNoise(cellrange,:)=NaN;
         
        drift_theta_w(cellrange,:)=NaN;
        drift_theta(cellrange,:)=NaN;
        
        driftspont(cellrange,:) = NaN;
              
        driftwpref(cellrange,:) = NaN;
        driftwbw(cellrange,:) = NaN ;
        
        driftF1F0(cellrange,:) = NaN;
        driftF0(cellrange,:) = NaN;
        %       driftorientfreq_all(cellrange,:)=field2array(drift, 'orientfreq_all');
        cvDSI(cellrange,:) = NaN; %%also circular variance measure change to "cv_dsi and cv_osi" in new compile programs
        cvOSI(cellrange,:)=NaN; 
        end
        
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

clear m ind x y t_lag STA STA1 STA_post STA1_post c crf crf1

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
    crf = wn(c).crf;
    crf1{c} = crf;

% plot
    end
    CRF(cellrange)=crf1;
    wn_spont(cellrange,:)= field2array_false(wn,'spont');
    
else
    STA_peak(cellrange)=NaN
%     crf1_post(cellrange)=NaN
    wn_spont(cellrange)=NaN
end
%close all
clear w

if exist ('wn_post','var')
    for w = 1:length(wn_post)
    STA_post = wn_post(w).sta;
    
    %%%Dtermine time point with maximial response
    [m ind] = max(abs(STA_post(:)-127));
    [x y t_lag] = ind2sub(size(STA_post),ind);
    
    STA1_post{w} = STA_post(:,:,t_lag)-128;

% figure
% imagesc(STA1{1,w}',[-64 64]); axis equal
    end

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

figure
plot(alldata(find(inh),5),alldata(find(inh),6),'ro');
hold on
plot(alldata(find(exc),5),alldata(find(exc),6),'ko');
hold on
plot(alldata(find(midnarrow),5),alldata(find(midnarrow),6),'bo');



for i = 1:size(driftA1,1)
    for j=1:size(driftA1,2)
        driftA1(i,j);
        driftA2(i,j);
        driftB(i,j);
        drift_theta_w(i,j);
        [OSI(i,j) DSI(i,j) width(i,j) peak(i,j)] = calculate_tuning(driftA1(i,j),driftA2(i,j),driftB(i,j),drift_theta_w(i,j));
        
        
    end
end
lyr = alldata(:,26)

%%% define responsive
resp = peak(:,2)>=2;

%%% define groups
DOI = GT'==3;
Lis=GT'==2;
EtOH = GT'==1;

%%%plot waveforms and STAs for pre drug treatment
use = find(resp & DOI);
for i=1:ceil(length(use)/8);
       figure
       for j= 1:8
 
    subplot(4,4,2*(j-1)+1);
    plot(wvform(use((i-1)*8+j),:));axis off

    subplot(4,4,2*(j-1)+2); 
    
    %plot(plotrange,psth_pinp(use((i-1)*8+j),plotrange));ylim([0 50]); xlim([min(plotrange) max(plotrange)]);
       
    if ~isempty(STA_peak{1,use((i-1)*8+j)})
    colormap jet
    imagesc(STA_peak{1,use((i-1)*8+j)}',[-64 64]);
    end
    title('pre DOI');
    %   title(sprintf('%0.0f',cells(use((i-1)*8+j))));
  set(gca,'Ytick',[]); set(gca,'Xtick',[])
       end
end


use = find( resp & Lis);
for i=1:ceil(length(use)/8);
       figure
       for j= 1:8
 
    subplot(4,4,2*(j-1)+1);
    plot(wvform(use((i-1)*8+j),:));axis off

    subplot(4,4,2*(j-1)+2); 
    
    %plot(plotrange,psth_pinp(use((i-1)*8+j),plotrange));ylim([0 50]); xlim([min(plotrange) max(plotrange)]);
       
    if ~isempty(STA_peak{1,use((i-1)*8+j)})
    colormap jet
    imagesc(STA_peak{1,use((i-1)*8+j)}',[-64 64]);
    end
    title('pre Lis');
    %   title(sprintf('%0.0f',cells(use((i-1)*8+j))));
  set(gca,'Ytick',[]); set(gca,'Xtick',[])
       end
end


%%%plot waveforms and STAs post drug treatment

use = find(resp & DOI);
for i=1:ceil(length(use)/8);
       figure
       for j= 1:8
 
    subplot(4,4,2*(j-1)+1);
    plot(wvform(use((i-1)*8+j),:));axis off

    subplot(4,4,2*(j-1)+2); 
    
    %plot(plotrange,psth_pinp(use((i-1)*8+j),plotrange));ylim([0 50]); xlim([min(plotrange) max(plotrange)]);
       
    if ~isempty(STA_peak_post{1,use((i-1)*8+j)})
    colormap jet
    imagesc(STA_peak_post{1,use((i-1)*8+j)}',[-64 64]);
    end
    title('post DOI');
    %   title(sprintf('%0.0f',cells(use((i-1)*8+j))));
  set(gca,'Ytick',[]); set(gca,'Xtick',[])
       end
end

use = find( resp & Lis);
for i=1:ceil(length(use)/8);
       figure
       for j= 1:8
 
    subplot(4,4,2*(j-1)+1);
    plot(wvform(use((i-1)*8+j),:));axis off

    subplot(4,4,2*(j-1)+2); 
    
    %plot(plotrange,psth_pinp(use((i-1)*8+j),plotrange));ylim([0 50]); xlim([min(plotrange) max(plotrange)]);
       
    if ~isempty(STA_peak_post{1,use((i-1)*8+j)})
    colormap jet
    imagesc(STA_peak_post{1,use((i-1)*8+j)}',[-64 64]);
    end
    title('post Lis');
    %   title(sprintf('%0.0f',cells(use((i-1)*8+j))));
  set(gca,'Ytick',[]); set(gca,'Xtick',[])
       end
end


resp = peak(:,2)>=2 & wn_spont1>0 ;

PlotSpontData(wn_spont1,wn_spont1_post,DOI,lyr,inh,resp,'DOI')
PlotSpontData(wn_spont1,wn_spont1_post,Lis,lyr,inh,resp,'Lis')


%%wt gain modulation of neurons
peak_run_p=nanmedian(peak(DOI & lyr<=4 & resp & ~inh,2))
peak_stat_p=nanmedian(peak(DOI & lyr<=4 & resp & ~inh,1))
gain_ind_p=(peak_run_p-peak_stat_p)/(peak_run_p+peak_stat_p)
sprintf('gain DOI %f',mean(gain_ind_p))

layerAgePlot_ratio_jlh(peak(:,1),peak(:,2),GT,lyr,inh,resp,{'run','stat'},'run vs stat');

%%% plot grating respons data
layerTreatmentPlot(peak(:,1),GT,lyr,inh,resp,'peak')

%%% waveforms of various cells
figure
plot(wvform(exc,:)','b');
hold on
plot(wvform(inh,:)','r'); hold on
plot(wvform(DOI,:)','g')


%generate STA for cells meeting criterion defined by l

l=find (~inh & DOI & resp);

for j=1:length(STA_peak)

    if ismember(j,l)
figure
if ~isempty(STA_peak{1,j})
    colormap jet
    imagesc(STA_peak{1,j}',[-64 64]);
end
    end
end

%%%spit out measures for cells of interest
driftF1F0(629,1)
OSI(629,1)
DSI(629,1)
driftwpref(629,1)










