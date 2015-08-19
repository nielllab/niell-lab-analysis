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

        
        afiles = {'NR5A_Pinping\12_22_14_pos_ctl\analysis_12_22_14_1st_clust.mat',...
           'NR5A_Pinping\2_25_15\analysis_2_25_15_rec1.mat',...
           'NR5A_Pinping\3_7_15\analysis_3_7_15.mat',...
           'NR5A_Pinping\3_11_15_wt\full\analysis_3_11_15_full.mat',...
           'NR5A_Pinping\3_13_15\analysis_3_15_15.mat',...
           'NR5A_Pinping\3_24_15\full\analysis_3_24_15.mat',...
           'NR5A_Pinping\4_9_15\Full\analysis_full.mat',...
           'NR5A_Pinping\4_10_15\analysis.mat',...
           'NR5A_Pinping\4_13_15\analysis.mat'};

  
   
    
    
    
    for i = 1:length(afiles)
        
        sessionNum=sessionNum+1;
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
        clear psth
        clear psth_power1
        clear psth_power2
        clear psth_power3
        
        
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
     
        
       
%                    
%                     
        
      
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
    
    STA_peak(cellrange)=STA1;
   
else
    STA_peak(cellrange)=NaN
end
%close all
%size(wvform)
        wvform(cellrange,:) = wv';
       
        %% get peristimulus histograms for PINPed units
       if exist('psth_power1', 'var'); 
        psth_pinp_1(cellrange,:)=psth_power1;
       else
           psth_pinp_1(cellrange,:)=NaN;
           
       end
       
       if exist('psth_power2', 'var'); 
        psth_pinp_2(cellrange,:)=psth_power2;
       else
           psth_pinp_2(cellrange,:)=NaN; 
       end
       
       if exist('psth_power3', 'var'); 
        psth_pinp_3(cellrange,:)=psth_power3;
       else
           psth_pinp_3(cellrange,:)=NaN;
           
       end
       
       if exist('psth'); 
        psth_pinp(cellrange,:)=psth;
       else
        psth_pinp(cellrange,:)=NaN;
       end
       
        histbins=histbins
    
        %get firing rate at all measured orients and SF, put into an array:
        %12 rows(orientations) by 7 columns(SpatialFreqs) for each cell
        
%         drift_Ori_Sf(cellrange,:) = arrayfun(@(x)(getfield(x,'orientfreq_all')),drift,'UniformOutput',false);
%        % drift_all(cellrange,:)=drift';
        
        if exist('locomotion');
           Vel(i,dataset)=arrayfun(@(x)(getfield(x,'mouseV')),locomotion,'UniformOutput',false);  
        end
       
    end %%% loop over adult vs EO


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
midnarrow = alldata(:,6)>0 & alldata(:,6)<4 & alldata(:,5)<12 &  alldata(:,5)>8.5;  %%% could analyze these specifically at some point
exc= alldata(:,5)>10 & alldata(:,6)>0 & alldata(:,6)<4;

figure
plot(alldata(find(inh),5),alldata(find(inh),6),'ro');
hold on
plot(alldata(find(exc),5),alldata(find(exc),6),'ko');
hold on

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
resp = peak(:,1)>=2

figure
plot(psth_pinp(lyr==4 & exc,:)')
set(gca,'xlim',[50 60])

psth_pinp_1=psth_pinp
psth_pinp_1(psth_pinp_1(:,:)>180)=0;

figure
plot(psth_pinp(lyr<=4 ,:)')



figure
plot(52:80,mean(psth_pinp(lyr==2&exc &resp ,52:80)),'m');
hold on
plot(52:80,mean(psth_pinp(lyr==3&exc & resp,52:80)),'b');
plot(52:80,mean(psth_pinp(lyr==4&exc & resp,52:80)),'g');
plot(52:80,mean(psth_pinp(lyr==5&exc & resp,52:80)),'k');
plot(52:80,mean(psth_pinp(inh & resp,52:80)),'r');
legend({'2','3','4','5','inh'})

%%% calculate baseline, evoked, and artifact by choosing windows
baseline = mean(psth_pinp(:,1:45),2);
baseStd = std(psth_pinp(:,1:45),[],2);
artifact = psth_pinp(:,50:52);
ev = mean(psth_pinp(:,53:54),2)
evoked1 = ev-baseline;

 
% ev = mean(psth_pinp(:,52:55),2) - mean(psth_pinp(:,[51 56 57]),2);
% evoked1 = ev 

% ev = max(psth_pinp(:,53:54),[],2)- min(psth_pinp(:,53:54),[],2);
% evoked = ev;

zscore1 =evoked1./baseStd;
zscore1(zscore1>10)=10;

%%% plot zscore vs evoked
figure
plot(zscore1,evoked1,'mo')
hold on
plot(zscore1(lyr==3),evoked1(lyr==3),'bo')
plot(zscore1(lyr==4),evoked1(lyr==4),'go')
plot(zscore1(inh),evoked1(inh),'ro');
xlabel('zscore'); ylabel('evoked')
legend({'all','lyr 3','lyr4','inh'})

%%% choose pinped
pinped = (zscore1>8& evoked1>8 ); 

figure
plot(45:80,psth_pinp(pinped ,45:80),'k');hold on
plot(45:80,psth_pinp(pinped & inh ,45:80),'r');hold on

figure
plot(wvform(exc & pinped & session'==5,:)','color','k');hold on
plot(wvform(inh & pinped & session'==5,:)','color','r');hold on

figure
plot(psth_pinp( pinped & exc & session'==5,:)')
xlim([50 60])


figure
bar([sum(exc & lyr==2 & pinped)/sum(exc & lyr==2) sum(exc & lyr==3 & pinped)/sum(exc & lyr==3) ...
    sum(exc & lyr==4&pinped)/sum(exc & lyr==4) sum(exc & lyr==5&pinped)/sum(exc & lyr==5) sum(inh&pinped)/sum(inh)]);
set(gca,'xticklabel',{'2','3','4','5','inh'},'ylim',[0 0.15])
title(sprintf('%d pinped',sum(pinped)))
ylabel('fraction pinped')

figure
plot(pinped(lyr==4),'*')
hold on
plot(GT(lyr==4)/3,'g')
set(gca,'ylim',[0.2 1.25])
xlabel('cell #')
legend({'lyr 4 pinped','genotype'})

figure
plot(psth_pinp(pinped,:)')
set(gca,'xlim',[50 60])

figure
plot(psth_pinp(~pinped & lyr==4,:)')
set(gca,'xlim',[50 60])

figure
plot(baseline,ev,'ko')
hold on
plot(baseline(pinped &lyr==4),ev(pinped&lyr==4),'go')
plot(baseline(pinped &lyr==3),ev(pinped&lyr==3),'bo')
plot(baseline(pinped&inh),ev(pinped&inh),'ro');
%axes([0 50  0 200])
legend('non','pinp lyr 4','pinp lyr3','pinp lyr2','pinp lyr5', 'pinp inh')
plot([0 50],[0 200])
xlabel('baseline'); ylabel('laser')
sprintf('%d pinped neurons total',sum(pinped))




figure
plot(baseline,ev,'ko')
axis equal; hold on
plot(baseline(pinped &GT'==1),ev(pinped&GT'==1),'ro');hold on
plot(baseline(pinped &GT'==2),ev(pinped&GT'==2),'bo');hold on
plot(baseline(pinped &GT'==3),ev(pinped&GT'==3),'go');hold on
legend('non','wt','N2A','N2B')
plot([0 40],[0 40])
xlabel('baseline'); ylabel('laser')


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
frac_resp(3) = sum(resp& lyr==4 & GT'==2 & pinped)/sum( lyr==4 & GT'==2  & pinped);

figure
bar(frac_resp); ylim([0 2]); ylabel('fraction resp >2'); title('lyr 4')
set(gca,'xticklabel',{'2B pinp','wt pinp','2A pinp'})

%%% fraction responsive
frac_resp(1) = sum(resp& lyr==4 & GT'==1 & ~pinped)/sum( lyr==4 & GT'==1 & ~pinped);
frac_resp(2) = sum(resp& lyr==4 & GT'==3  & ~pinped)/sum( lyr==4 & GT'==3  & ~pinped);
frac_resp(3) = sum(resp& lyr==4 & GT'==2 & ~pinped)/sum( lyr==4 & GT'==2  & ~pinped);

figure
bar(frac_resp); ylim([0 2]); ylabel('fraction resp >2'); title('lyr 4')
set(gca,'xticklabel',{'2B non pinp','wt non pinp','2A pinp'})


%%wt gain modulation of pinped neurons
peak_run_p=nanmedianMW(peak(wt & lyr<=4 & pinped & resp & exc,2))
peak_stat_p=nanmedianMW(peak(wt & lyr<=4 & pinped & resp & exc,1))
gain_ind_p=(peak_run_p-peak_stat_p)/(peak_run_p+peak_stat_p)

sprintf('gain wt pinped %f',mean(gain_ind_p))

peak_run=nanmedianMW(peak(wt & lyr ==4 & ~pinped & resp & exc,2))
peak_stat=nanmedianMW(peak(wt &lyr==4& ~pinped & resp & exc,1))
gain_ind=(peak_run-peak_stat)/(peak_run+peak_stat)
sprintf('gain wt not pinped %f',mean(gain_ind))

figure
bar([peak_stat_p peak_run_p; peak_stat peak_run])
title('peak resp wt')
legend('stationary','run');
set(gca,'xticklabel',{'pinped','non-pinped'})


%%NR2B gain modualtion of pinped neurons

peak_run_p=nanmedianMW(peak(N2B & lyr<=4 & pinped & resp & exc,2))
peak_stat_p=nanmedianMW(peak(N2B & lyr<=4 & pinped & resp & exc,1))
gain_ind_p=(peak_run_p-peak_stat_p)/(peak_run_p+peak_stat_p)

sprintf('gain n2b pinped %f',mean(gain_ind_p))

peak_run=nanmedianMW(peak(N2B & lyr==4 & ~pinped & resp & exc,2))
peak_stat=nanmedianMW(peak(N2B & lyr==4 & ~pinped & resp & exc,1))
gain_ind=(peak_run-peak_stat)/(peak_run+peak_stat)
sprintf('gain n2b not pinped %f',mean(gain_ind))
figure
bar([peak_stat_p peak_run_p; peak_stat peak_run])
title('peak resp NR2B KO')
legend('stationary','run');
set(gca,'xticklabel',{'pinped','non-pinped'})


%layerAgePlot_ratio_jlh(peak(:,1),peak(:,2),GT,lyr,inh,resp,{'run','stat'},'run vs stat');



%%% plot grating respons data

lyr34 = lyr==3 | lyr==4;

plotPinpData(driftspont,wt,lyr==4,pinped,inh,1)
plot([0 20],[0 20])
title('spont')

barPinp(driftspont,wt,N2B,N2A,lyr==4,pinped,inh,1)
ylabel('spont')

plotPinpData(peak,wt,lyr==4,pinped,inh,1)
plot([0 20],[0 20])
title('peak')

barPinp(peak(:,2),wt,N2B,N2A, lyr==4,pinped,inh,1)
ylabel('peak')



plotPinpData(driftF1F0,wt,lyr,pinped,inh,resp)
plot([0 2],[0 2])
title('F1F0')

barPinp(driftF1F0(:,1),wt,N2A,N2B,lyr==4,pinped,inh,resp)
ylabel('F1F0')


plotPinpData(cvOSI,wt,lyr==4,pinped,inh,resp)
plot([0 1],[0 1])
title('cvOSI')

barPinp(cvOSI(:,1),wt,N2B,N2A,lyr==4,pinped,inh,resp)
ylabel('cvOSI')

plotPinpData(OSI,wt,lyr,pinped,inh,resp)
plot([0 1],[0 1])
title('OSI')

barPinp(OSI(:,1),wt,N2A,N2B,lyr==4,pinped,inh,resp)
ylabel('OSI')

plotPinpData(cvDSI,wt,lyr==4,pinped,inh,resp)
plot([0 1],[0 1])
title('cv DSI')

barPinp(cvDSI(:,1),wt,N2B,N2A,lyr==4,pinped,inh,resp)
ylabel('cvDSI')

barPinp(DSI(:,1),wt,N2B,N2A,lyr==4,pinped,inh,resp)
ylabel('cvDSI')


driftwpref(driftwpref==0) = 0.005; %tranform SF pref = 0 to 0.005 in order to plot on log scale
SF_pref_stat = resp & driftwpref(:,1) >=0.005 & driftwpref(:,1)<=0.40;

barPinp(driftwpref(:,1),wt,N2B,N2A,lyr==4,pinped,inh,SF_pref_stat)
ylabel('SFpref')

drift_theta_1=size(drift_theta);
drift_theta_1=(drift_theta*180)/pi;
drift_theta_1(drift_theta_1>330)=0;
%%preferred oprientation top 30prct
layerAgePlot_pref_Orient(drift_theta_1(:,1),GT,lyr,inh,resp,{'pref Orient Cardinal' 'Prct total'},'Prefered Orientation Stationary');

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
k=find(pinped &exc & N2A)
l=find(~pinped& exc  & wt)
n=find(pinped&inh)
% for j=length(I)
% figure
% imagesc(STA_peak{1,12});
% end

h=ismember(k,I)

%generate STA for pinped cells

for j=1:length(STA_peak)

    if ismember(j,k)
figure
if ~isempty(STA_peak{1,j})
    colormap jet
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







