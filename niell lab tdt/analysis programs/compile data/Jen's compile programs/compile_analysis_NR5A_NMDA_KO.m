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


for dataset = 1:3  %%% control ("wt") NR5A1-cre/CHR2 animals vs. NR2A deleted NR5A1-cre
    
    if dataset ==1
        
        
        afiles = { 'NR5A_Pinping\1_13_15_pos_ctl\1_13_15_analysis_pos_ctl.mat',...
            'NR5A_Pinping\12_22_14_pos_ctl\analysis_12_22_14_1st_clust.mat'};

    elseif dataset==2
%         
        afiles = {'NR2A\KO\8_12_14_NR2A_KO\Rec2\analysis_8_12_14_rec2.mat',...
            'NR2A\KO\9_19_14\rec2\analysis_9_19_14_NR2A_rec2.mat',...
            'NR2A\KO\10_8_14\full clustering_rec1\analysis_10_8_14_NR2A_KO_rec1.mat'};
%         
    elseif dataset==3
        afiles = {'NR2B\1_24_15_KO\analysis_1_24_15_KO.mat',...
            'NR2B\12_8_14_homo\analysis_12_8_14_NR2BKO.mat',...
            'NR2B\1_26_15_NR2B_KO_NR5A\analysis_1_26_15_NR2B_KO_NR5A1.mat',...
            'NR2B\2_2_15\analysis_2_2_15_NR2B_KO.mat'};
        
    end
    
    
    
    for i = 1:length(afiles)
        
        clear params
        clear wn 
        clear wn_movement
        clear LFP_movement
        clear bars
        clear wave_all
        clear rf_width
        clear locomotion
        clear drift
        clear layer
        
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
%          for c=1:length(cells);
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
% %             
% % %             
%         end
%         
%         SNRall(cellrange)=snr;
%         meanWavesAll(cellrange,:,:) = meanwaves;
%         stdWavesAll(cellrange,:,:) = stdwaves;
        
        if exist('params','var');

            
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
        
%         if exist ('wn_movement','var')
%            
%             stopCRF(cellrange)=field2array_false(wn_movement,'stopCRF');
%             moveCRF(cellrange)=field2array_false(wn_movement,'moveCRF');
%         else
%             stopCRF(cellrange)=NaN;
%             moveCRF(cellrange)=NaN;
%         end

clear m ind x y t_lag STA1

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
   
else
    STA_peak(cellrange)=NaN
end
%close all
%size(wvform)
        wvform(cellrange,:) = wv';
       
        %% get peristimulus histograms for PINPed units
        psth_pinp(cellrange,:)=psth;
        histbins=histbins
     
        %get firing rate at all measured orients and SF, put into an array:
        %12 rows(orientations) by 7 columns(SpatialFreqs) for each cell
        
%         drift_Ori_Sf(cellrange,:) = arrayfun(@(x)(getfield(x,'orientfreq_all')),drift,'UniformOutput',false);
%        % drift_all(cellrange,:)=drift';
        
        if exist('locomotion');
           Vel(i,dataset)=arrayfun(@(x)(getfield(x,'mouseV')),locomotion,'UniformOutput',false);  
        end
       
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


inh = alldata(:,6)>0 & alldata(:,6)<1.6  & alldata(:,5)<8.0;  %%% directly based on wvform; k-means includes another inh group?
midnarrow = alldata(:,6)>0 & alldata(:,6)<4 & alldata(:,5)<12 &  alldata(:,5)>8.5;  %%% could analyze these specifically at some point
exc= alldata(:,5)>10 & alldata(:,6)>0 & alldata(:,6)<4;

figure
plot(alldata(find(inh),5),alldata(find(inh),6),'ro');
hold on
plot(alldata(find(exc),5),alldata(find(exc),6),'ko');
hold on
%plot(alldata(find(pinp),5),alldata(find(pinp),6),'go');

figure
plot(wvform(exc,:)','color','k');hold on
plot(wvform(inh,:)','color','r');hold on


for i = 1:size(driftA1,1)
    for j=1:size(driftA1,2)
        driftA1(i,j);
        driftA2(i,j);
        driftB(i,j);
        drift_theta_w(i,j);
        [OSI(i,j) DSI(i,j) width(i,j) peak(i,j)] = calculate_tuning(driftA1(i,j),driftA2(i,j),driftB(i,j),drift_theta_w(i,j));
        
        
    end
end

%%% calculate baseline, evoked, and artifact by choosing windows
baseline = mean(psth_pinp(:,5:45),2);
baseStd = std(psth_pinp(:,5:45),[],2);
artifact = psth_pinp(:,50);
ev = mean(psth_pinp(:,51:55),2);
evoked = ev-baseline;

zscore =evoked./baseStd;
zscore(zscore>10)=10;
lyr=alldata(:,26);

%%% plot zscore vs evoked
figure
plot(zscore,evoked,'mo')
hold on
plot(zscore(lyr==3),evoked(lyr==3),'bo')
plot(zscore(lyr==4),evoked(lyr==4),'go')
plot(zscore(lyr==5),evoked(lyr==5),'ro');
plot(zscore(lyr==6),evoked(lyr==6),'yo');
legend('lyr 2','lyr3','4','5','6')
xlabel('zscore'); ylabel('evoked')

% artscore = (artifact - baseline)./baseStd;
% artscore(artscore>10)=10;
% 
% figure
% hist(artscore',0:10)

%%% choose pinped
pinped = (zscore>2& evoked>5); 

figure
plot(baseline,ev,'ko')
axis equal;hold on
plot(baseline(pinped &lyr==4),ev(pinped&lyr==4),'go');hold on
plot(baseline(pinped &lyr==3),ev(pinped&lyr==3),'bo');hold on
plot(baseline(pinped &lyr<=2),ev(pinped&lyr<=2),'mo');hold on
plot(baseline(pinped &lyr==5),ev(pinped&lyr==5),'ro');hold on
plot(baseline(pinped&inh),ev(pinped&inh),'rs');

legend('non','pinp lyr 4','pinp lyr3','pinp lyr2','pinp lyr5', 'pinp inh')
plot([0 40],[0 40])
xlabel('baseline'); ylabel('laser')
sprintf('%d pinped neurons total',sum(pinped))


figure
plot(baseline,ev,'ko')
axis equal; hold on
plot(baseline(pinped &GT'==1),ev(pinped&GT'==1),'ro');hold on
plot(baseline(pinped &GT'==2),ev(pinped&GT'==2),'go');hold on
plot(baseline(pinped &GT'==3),ev(pinped&GT'==3),'bo');hold on
legend('non','pinp 2B KO','pinp 2A KO','pinp wt')
plot([0 40],[0 40])
xlabel('baseline'); ylabel('laser')



%%% define responsive
resp = driftpeak(:,2)>2 & driftpeak(:,1)>2;

%%% # spikes evoked
nbins=8;
ev_spikes = mean(psth_pinp(pinped,54 + (1:nbins)),2)-baseline(pinped);
ev_spikes = nbins*ev_spikes/1000;
figure
hist(ev_spikes)
xlabel('# evoked spikes - pinped');

nbins=8;
ev_spikes = mean(psth_pinp(pinped,54 + (1:nbins)),2)-baseline(pinped);
ev_spikes = nbins*ev_spikes/1000;
h1 = hist(ev_spikes,-0.05:0.02:0.25)/length(ev_spikes);
ev_spikes = mean(~psth_pinp(~pinped,54 + (1:nbins)),2)-baseline(~pinped);
ev_spikes = nbins*ev_spikes/1000;
h2 = hist(ev_spikes,-0.05:0.02:0.25)/length(ev_spikes);
figure
bar((-0.05:0.02:0.25)+0.01,[h1; h2]')

%%% define wt
wt = GT'==3;
N2A=GT'==2;
N2B = GT'==1;

%%% fraction responsive
frac_resp(1) = sum(resp& lyr==4 & GT'==1 & pinped)/sum( lyr==4 & GT'==1 & pinped);
frac_resp(2) = sum(resp& lyr==4 & GT'==3  & pinped)/sum( lyr==4 & GT'==3  & pinped);
figure
bar(frac_resp); ylim([0 2]); ylabel('fraction resp >2'); title('lyr 4')
set(gca,'xticklabel',{'2B pinp','wt pinp'})


%%wt gain modulation of pinped neurons
peak_run_p=nanmedianMW(peak(wt & lyr<=4 & pinped & resp & exc,2))
peak_stat_p=nanmedianMW(peak(wt & lyr<=4 & pinped & resp & exc,1))
gain_ind_p=(peak_run_p-peak_stat_p)/(peak_run_p+peak_stat_p)


peak_run=nanmedianMW(peak(wt & lyr <=4 & ~pinped & resp & exc,2))
peak_stat=nanmedianMW(peak(wt &lyr<=4& ~pinped & resp & exc,1))
gain_ind=(peak_run-peak_stat)/(peak_run+peak_stat)


%%NR2B gain modualtion of pinped neurons

peak_run_p=nanmedianMW(peak(N2B & lyr<=4 & pinped & resp & exc,2))
peak_stat_p=nanmedianMW(peak(N2B & lyr<=4 & pinped & resp & exc,1))
gain_ind_p=(peak_run_p-peak_stat_p)/(peak_run_p+peak_stat_p)


peak_run=nanmedianMW(peak(N2B & lyr<=4 & ~pinped & resp & exc,2))
peak_stat=nanmedianMW(peak(N2B & lyr<=4 & ~pinped & resp & exc,1))
gain_ind=(peak_run-peak_stat)/(peak_run+peak_stat)

layerAgePlot_ratio_jlh(peak(:,1),peak(:,2),GT,lyr,inh,resp,{'run','stat'},'run vs stat');



%%% plot grating respons data

plotPinpData(driftspont,wt,1,pinped,inh,1)
plot([0 20],[0 20])
title('spont')

plotPinpData(peak,wt,1,pinped,inh,1)
plot([0 20],[0 20])
title('peak')

plotPinpData(driftF1F0,wt,1,pinped,inh,resp)
plot([0 2],[0 2])
title('F1F0')

plotPinpData(cvOSI,N2B,1,pinped,inh,resp)
plot([0 1],[0 1])
title('cvOSI')

plotPinpData(OSI,N2B,1,pinped,inh,resp)
plot([0 1],[0 1])
title('OSI')



plotPinpData(driftF1F0,N2B,1,pinped,inh,resp)
plot([0 2],[0 2])
title('F1F0 2B pinp')

plotPinpData(DSI,N2B,1,pinped,inh,resp)
plot([0 1],[0 1])
title('DSI')


plotPinpData(cvDSI,wt,lyr<=4,pinped,inh,resp)
plot([0 1],[0 1])
title('cv DSI')

plotPinpData(driftpeak,GT'==1,lyr==4,pinped,inh,resp)
hold on
plot([0 50],[0 50])
title('peak')

%%% show psth of all pinped cells
%%% lots of figures but useful

% pinplist = find(pinped);
% for i=1:length(pinplist)
% c = pinplist(i);
% 
%     figure
% plot(psth_pinp(c,:)')
% hold on
% plot([0 100],[baseline(c) baseline(c)],'g')
% plot([0 100],[baseline(c)+baseStd(c)*3 baseline(c)+baseStd(c)*3],'r')
% title(sprintf('evoked = %f  zscore = %f layer = %d inh = %d gt = %d',evoked(c),zscore(c),lyr(c),inh(c),GT(c)))
% end


%%% waveform of pinped vs non
figure
plot(wvform(exc,:)','b');
hold on
plot(wvform(inh,:)','r'); hold on
plot(wvform(pinped,:)','g')
figure
plot(wvform(pinped,:)','g')
title('pinped')

I=find(inh)
k=find(pinped & wt)
l=find(~pinped&N2B)
n=find(pinped&inh)
% for j=length(I)
% figure
% imagesc(STA_peak{1,12});
% end

h=ismember(k,I)

%generate STA for pinped cells
for j=1:length(STA_peak)

    if ismember(j,l)
figure
if ~isempty(STA_peak{1,j})
    imagesc(STA_peak{1,j},[-64 64]);
end
    end
end

%generate STA for inhibitory cells
for j=1:length(STA_peak)

    if ismember(j,I)
figure
imagesc(STA_peak{1,j},[-64 64]);
    end
end

figure
plot(wvform(inh,:)','color','r');hold on
for j=1:length(STA_peak)

    if ismember(j,k)
plot(wvform(j,:)','color','g'); hold on
    end
end




figure
subplot(1,2,1);
pie([sum(wt&pinped&~inh) sum(wt&pinped&inh)],{'broad','narrow'})
title('wt pinped')


figure
subplot(1,2,1);
pie([sum(N2B&pinped&~inh) sum(N2B&pinped&inh)],{'broad','narrow'})
title('N2B pinped')

figure
subplot(1,2,1);
pie([sum(N2A&pinped&~inh) sum(N2A&pinped&inh)],{'broad','narrow'})
title('N2A pinped')


% subplot(1,2,2);
% pie([sum(~wt&pinped&~inh) sum(~wt&pinped&inh)],{'broad','narrow'})
% title('2A/2B ko pinped')

%%%                         Receptive field STA data

%%%Nx data
clear f x f2 x1
good_STA =STA_exp_var>=0.6 & STA_ny<0.39;
responsive_stat = peak(:,1)>=2;  % firing rate (responsiveness) criteria for whether cells enter subsequent statistical analysis







