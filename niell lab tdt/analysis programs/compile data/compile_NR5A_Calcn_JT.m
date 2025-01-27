close all
clear all
BatchEphys_Johanna

use =  find( ~cellfun(@isempty,{files.prefPinp})  )
%use= use(1); %%% just analyze first one for now!

%use = find(strcmp({files.expt},'1_14_16'))


n=0;
for i = 1:length(use)
    afile = [pathname files(use(i)).path  files(use(i)).analysisfile ];
    clustfile = [pathname files(use(i)).path files(use(i)).clusterfile] ;
    
    %%% get cell type based on waveform
    [inh mid] = getWaveform(clustfile,afile,0);
    nc = length(inh); cellrange = n+1:n+nc;
    inhAll(cellrange) = inh;
    
    %%% get layer info
    redo=0;
    layer = getLayer(clustfile,afile,files(use(i)).tip_loc_1,files(use(i)).tip_loc_2,files(use(i)).angle,redo);
    allLayer(cellrange) = layer.units;
  
    %%% experiment parameters  
    condition(cellrange)=files(use(i)).condition;
    age(cellrange) = files(use(i)).age;
    sex(cellrange) = files(use(i)).sex;
    session(cellrange) = i;
    
    %%% get pinping info
    trainPeriod = 5; freq = 5; showImg=0; makePDF=1; redo=0;
    pinp = getPinp(clustfile,afile,files(use(i)).prefPinp,trainPeriod,freq,showImg,makePDF,redo);
    pinped(cellrange) = pinp.pinped; %% append
       
    %%% get drifting grating response
    redo = 0;
    drift = getDrift_mv(clustfile,afile,files(use(i)).blockDrift,redo);
    drift_orient(cellrange,:,:)=drift.orient_tune;
    drift_sf(cellrange,:,:) = drift.sf_tune;
    drift_spont(cellrange,:) = drift.interSpont;
    drift_osi(cellrange,:) = drift.cv_osi;
    drift_F1F0(cellrange) = drift.F1F0;   
    
     %%% get white noise response
%     redo = 1;
%     periodFrames=10;%number of frames with increasing or decreasing contrast,
%     wn_mv= getWn_mv(clustfile,afile,files(use(i)).blockWn,redo,periodFrames);
   
    
    %%% get grating speed
        spd = getSpeed(clustfile,afile,files(use(i)).blockDrift,0);
        speedHistDrift(i,:) = hist(spd.v,0.5:1:100)/length(spd.v);
        speedTrace{i}=spd.v;
    
    n= n+nc;
end
inh = inhAll;
exc = ~inh;

pinpedI=(pinped & inh);
pinped = (pinped & exc); %%% pinped units should be excitatory

maxresp = max(squeeze(drift_orient(:,:,2)),[],2);  %%% max grating response moving
responsive = maxresp'>2;

%%% examples of analyses

%%% OSI by condition
osibins = 0.05:0.1:1;
figure
hist(drift_osi(condition==1 &responsive),osibins);
hold on
hist(drift_osi(condition==2 & responsive),osibins);
legend('condition 1','condition 2')

%%% OSI pinped and non-pinped
figure
hist(drift_osi(~pinped & condition==2 & exc &responsive),osibins);
hold on
hist(drift_osi(pinped & exc & responsive & condition==2),osibins);
legend('non-pinped','pinped')

%%% F1)1 pinped and non-pinped
f1f0bins = 0:0.1:2;
figure
hist(drift_F1F0(~pinped & exc &responsive),f1f0bins);
hold on
hist(drift_F1F0(pinped & responsive),f1f0bins);
legend('non-pinped','pinped')



%%% plot OSI tuning curves
%%% denotes pinped
for s = 1:max(session)

    figure
    set(gcf,'Name',sprintf('%s condition %d',files(use(s)).expt,files(use(s)).condition));
    useN = find(session ==s)
    for i = 1:length(useN)
        np = ceil(sqrt(length(useN)));
        subplot(np,np,i);
        hold on
        plot(drift_orient(useN(i),:,1),'Color',[1 0 0]);  plot(drift_orient(useN(i),:,2),'Color',[0 1 0]); %pre sal & DOI mv & stat
        plot([1 12], [1 1]*drift_spont(useN(i),1),':','Color',[1 0 0]);  plot([1 12], [1 1]*drift_spont(useN(i),2),':','Color',[0 1 0]);
        yl = get(gca,'Ylim'); ylim([0 max(yl(2),10)]); xlim([0.5 12.5])
        if pinped(useN(s))
            plot(6,8,'b*')
        end
    end

end

%%% speed distribution
figure
plot(speedHistDrift);
axis([0 40 0 1])
xlabel('speed'); ylabel('fraction')

figure
for i = 1:length(speedTrace);
    plot(0.1*(1:length(speedTrace{i}))/60,speedTrace{i});
end
xlabel('min'); ylabel('cm/sec')

keyboard





%%% old batch mode
%%% lots of useful figures!
%%% should be updated to new batch mode



N =0; cells=0;  all_img_STA={};PINPed=0; STA_peak={};stopCRF={}; moveCRF={};
sessionNum=0;

%takes the fields within your structure and transforms them into arrays so
%that you can access their information more readily

Cond=field2array(files,'condition');
Age=field2array(files,'age');
Sex=field2array(files,'sex');

idx_ctl=find(Cond==1);
idx_KO=find(Cond==2);


for dataset = 1:2  %%% control ("wt") NR5A1-cre/CHR2 animals vs. KO NR5A1-cre+ cells
    
    if dataset ==1
        
        for j=1:length(Cond);
            if  ismember(j,idx_ctl);
                afiles{j}=files(j).analysisfile;
                apath{j}=files(j).path;
                analysisfile{j}=[apath{j} afiles{j}];
                
            end
        end
        
    elseif dataset==2
        for j=1:length(Cond);
            if  ismember(j,idx_KO);
                afiles{j}=files(j).analysisfile;
                apath{j}=files(j).path;
                analysisfile{j}=[apath{j} afiles{j}];
                
            end
        end
        
    end
    
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
        
        %         h={files(1).path}
        %         l={files(1).analysisfile}
        %         m=[files(1).path files(1).analysisfile]
        %%
        load(analysisfile{i})
        
        clusterfilename
        
        afiles(i)
        
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
                save([apath{i} afiles{i}],'clusterfilename','-append');
            end
        end
        clusterFile
        try
            load(clusterFile,'wave_all');
        catch
            display('no cluster file')
        end
        
        %%
        
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
        
        
        GT(cellrange)=3-dataset;
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
        % % %
        % % % %
        %         end
        %
        %         SNRall(cellrange)=snr;
        %         meanWavesAll(cellrange,:,:) = meanwaves;
        %         stdWavesAll(cellrange,:,:) = stdwaves;
        
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
        
        
        
        clear m ind x y t_lag STA1 c crf crf1
        
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
                    %wn(c).spont = 0;
                end
                crf = wn(c).crf;
                crf1{c} = crf;
            end
            
            CRF(cellrange)=crf1;
        else
            STA_peak(cellrange)=NaN
            CRF(cellrange)=NaN
        end
        %close all
        %size(wvform)
        wvform(cellrange,:) = wv';
        
        %% get peristimulus histograms for PINPed units
        %        if exist('psth_power2','var');
        %         psth_pinp(cellrange,:)=psth_power2;
        %        else
        %            psth_pinp(cellrange,:)=psth;
        %        end
        
        psth_pinp(cellrange,:)=psth;
        histbins=histbins
        
        %get firing rate at all measured orients and SF, put into an array:
        %12 rows(orientations) by 7 columns(SpatialFreqs) for each cell
        
        %         drift_Ori_Sf(cellrange,:) = arrayfun(@(x)(getfield(x,'orientfreq_all')),drift,'UniformOutput',false);
        %        % drift_all(cellrange,:)=drift';
        
        if exist('locomotion');
            Vel(i,dataset)=arrayfun(@(x)(getfield(x,'mouseV')),locomotion,'UniformOutput',false);
        end
        
    end %%% loop over conditions
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

wave=alldata(:,7:25);

inh = alldata(:,6)>0 & alldata(:,6)<1.6  & alldata(:,5)<7.5  %%% directly based on wvform; k-means includes another inh group?
midnarrow = alldata(:,6)>0 & alldata(:,6)<4 & alldata(:,5)<10 &  alldata(:,5)>7.5;  %%% could analyze these specifically at some point
exc= alldata(:,5)>10 & alldata(:,6)>0 & alldata(:,6)<4;

lyr = alldata(:,26)

figure
plot(alldata(find(inh),5),alldata(find(inh),6),'ro');
hold on
plot(alldata(find(exc),5),alldata(find(exc),6),'ko');
hold on
plot(alldata(find(midnarrow),5),alldata(find(midnarrow),6),'bo');


%plot(alldata(find(pinp),5),alldata(find(pinp),6),'go');

%junk = (psth_pinp(:,52)>50);
junk= psth_pinp(:,52)>50;

figure
plot(psth_pinp(lyr<=4 & junk,:)')
set(gca,'xlim',[50 60])

figure
%plot(wvform(junk,:)','color','k');hold on
plot(wvform(inh &~junk ,:)','color','r');hold on
plot(wvform(midnarrow&~junk ,:)','color','b');hold on
plot(wvform(exc & ~junk ,:)','color','g');hold on
plot(wvform(junk,:)','color','k');hold on


for i = 1:size(driftA1,1)
    for j=1:size(driftA1,2)
        driftA1(i,j);
        driftA2(i,j);
        driftB(i,j);
        drift_theta_w(i,j);
        [OSI(i,j) DSI(i,j) width(i,j) peak(i,j)] = calculate_tuning(driftA1(i,j),driftA2(i,j),driftB(i,j),drift_theta_w(i,j));
        
        
    end
end
lyr = alldata(:,26);

%%% define responsive
resp = peak(:,1)>=2 %peak firing rate during gratings

wn_resp=0 % peak response to max contrast (max CRF at "10")
for i=1:length(CRF)
    if CRF{:,i}(10)>20;
        wn_resp(i)=1;
    else
        wn_resp(i)=0;
    end
end

plotrange = 50:80;
figure
plot(plotrange,mean(psth_pinp(lyr==2 &~inh & ~junk,plotrange)),'m');
hold on
plot(plotrange,mean(psth_pinp(lyr==3 &~inh &~junk,plotrange)),'b');
plot(plotrange,mean(psth_pinp(lyr==4 & ~inh &~junk,plotrange)),'g');
plot(plotrange,mean(psth_pinp(lyr==5 & ~inh & ~junk,plotrange)),'k');
plot(plotrange,nanmean(psth_pinp(inh & ~junk,plotrange)),'r');
legend({'2','3','4','5','inh'})


%%% calculate baseline, evoked, and artifact by choosing windows
baseline = mean(psth_pinp(:,5:45),2);
baseStd = std(psth_pinp(:,5:45),[],2);

%ev = mean(psth_pinp(:,53:55),2) - mean(psth_pinp(:,[52 56 57]),2);
ev = max(psth_pinp(:,53:54),[],2);
evoked = ev- baseline;

zscore =evoked./baseStd;
zscore(zscore>20)=20;

%%% plot zscore vs evoked
figure
plot(zscore(zscore>0),evoked(zscore>0),'mo')
hold on
plot(zscore(session'==20 & zscore>0),evoked(session'==20 & zscore>0),'bo')
plot(zscore(lyr==4& zscore>0),evoked(lyr==4& zscore>0),'go')
plot(zscore(lyr==5& zscore>0),evoked(lyr==5& zscore>0),'ko')
plot(zscore(inh& zscore>0),evoked(inh& zscore>0),'ro');
xlabel('zscore'); ylabel('evoked')
legend({'all','lyr 3','lyr4','lyr5','inh'})

%%% choose pinped
resp = peak(:,1)>=2;

pinped = (zscore>8& evoked>26  & ~inh  );
pinped_resp_grat = (zscore>8& evoked>26  & ~inh  & resp);
sum(pinped)
sum(pinped_resp_grat)
%%% define wt
wt = GT'==3;
N2A=GT'==2;
N2B = GT'==1;

frac_responsive_wt=(sum(pinped_resp_grat & wt)/sum(pinped & wt))
frac_responsive_N2B=(sum(pinped_resp_grat & N2B)/sum(pinped & N2B))
frac_responsive_N2A=(sum(pinped_resp_grat & N2A)/sum(pinped & N2A))

sum(pinped& N2A)
sum(pinped&N2B)
sum(pinped&wt)

sum(pinped_resp_grat& N2A)
sum(pinped_resp_grat&N2B)
sum(pinped_resp_grat&wt)
% for c = 1:length(CRF)
%     STA = wn(w).sta;
%
%     %%%Dtermine time point with maximial response
%     [m ind] = max(abs(STA(:)-127));
%     [x y t_lag] = ind2sub(size(STA),ind);
%
%     STA1{w} = STA(:,:,t_lag)-128;
%
%  end

use = find(wt & pinped);
clear crf
num_plot=ceil(length(use)/3)
for i = 1:ceil(length((use))/8)
    
    
    scrsz = get(0,'ScreenSize');
    figure('Position',[500 0 scrsz(3)/4 scrsz(4)])
    
    for j= 1:min(8,length(use)-(i-1)*8)
        
        if inh(use((i-1)*8+j))
            col = 'r';
        else col = 'b';
        end
        
        subplot(8,3,3*(j-1)+1);
        plot(wvform(use((i-1)*8+j),:),col);axis off
        title(['OS ',num2str(OSI(use((i-1)*8+j)))])
        subplot(8,3,3*(j-1)+2);
        if ~isempty(STA_peak{1,use((i-1)*8+j)})
            colormap jet
            imagesc(STA_peak{1,use((i-1)*8+j)}',[-64 64]);
            title(['SF pref ',num2str(driftwpref(use((i-1)*8+j)))])
        end
        
        
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

% figure
% imagesc(STA_peak{1,1263}',[-64 64]); axis equal
%     set(gca,'Ytick',[]); set(gca,'Xtick',[])

plotrange = 50:80;

for i=1:length(use)/8;
    figure
    for j= 1:8
        
        subplot(4,4,2*(j-1)+1);
        plot(wvform(use((i-1)*8+j),:));axis off
        
        subplot(4,4,2*(j-1)+2);
        plot(plotrange,psth_pinp(use((i-1)*8+j),plotrange));ylim([0 50]); xlim([min(plotrange) max(plotrange)]);
        hold on; plot([54 54], [ 0 50],'g')
        hold on; plot([52 52], [ 0 50],'r')
        title(sprintf('%0.0f %0.0f',evoked(use((i-1)*8+j)),zscore(use((i-1)*8+j))));
        set(gca,'Ytick',[]); set(gca,'Xtick',[])
    end
end


figure
plot(45:80,psth_pinp(pinped & ~junk,45:80),'g');hold on
plot(45:80,psth_pinp(inh ,45:80),'r');hold on
plot(45:80,psth_pinp(pinped & ~inh & ~junk ,45:80),'g');hold on
plot([50 50], [ 0 600],'g');hold on; plot([51 51], [ 0 600],'g')

% figure
% plot(psth_pinp( pinped & ~junk & ~inh ,:)','g');
% hold on
% plot(psth_pinp( pinped & ~junk & inh ,:)','r');
% xlim([50 60])

figure
bar([sum(~inh & lyr==2 & pinped_resp_grat & ~junk)/sum(~inh & lyr==2 & resp & ~junk) sum(~inh & lyr==3 & pinped_resp_grat & ~junk)/sum(~inh & lyr==3 & resp & ~junk) ...
    sum(~inh & lyr==4&pinped_resp_grat & ~junk)/sum(~inh & lyr==4 & resp & ~junk) sum(~inh & lyr==5& pinped_resp_grat & ~junk)/sum(~inh & lyr==5 & resp& ~junk) sum(inh & pinped_resp_grat & ~junk)/sum(inh& resp & ~junk)]);
set(gca,'xticklabel',{'2','3','4','5','inh'},'ylim',[0 0.60])
title(sprintf('%d pinped',sum(pinped_resp_grat)))
ylabel('fraction pinped')

figure
plot(pinped_resp_grat(lyr==4 & resp),'*')
hold on
plot(GT(lyr==4 & resp)/3,'g')
set(gca,'ylim',[0.2 1.25])
xlabel('cell #')
legend({'lyr 4 pinped','genotype'})


figure
plot(baseline(~junk),ev(~junk),'ko')
hold on
plot(baseline(pinped_resp_grat),ev(pinped_resp_grat),'go');hold on
plot(baseline(pinped_resp_grat &lyr==3),ev(pinped_resp_grat&lyr==3),'bo');hold on
plot(baseline(pinped_resp_grat &lyr==5),ev(pinped_resp_grat&lyr==5),'mo');hold on
plot(baseline(pinped_resp_grat&inh),ev(pinped_resp_grat&inh),'ro');
legend('non','pinp lyr 4','pinp lyr3','pinp lyr5', 'pinp inh')
plot([0 25],[0 25])
xlabel('baseline'); ylabel('laser')
sprintf('%d pinped neurons total',sum(pinped_resp_grat))

figure
plot(baseline(~junk),ev(~junk),'ko')
hold on
plot(baseline(pinped),ev(pinped),'go');hold on
plot(baseline(pinped &lyr==3),ev(pinped&lyr==3),'bo');hold on
plot(baseline(pinped &lyr==5),ev(pinped&lyr==5),'mo');hold on
plot(baseline(pinpedt&inh),ev(pinped&inh),'ro');
legend('non','pinp lyr 4','pinp lyr3','pinp lyr5', 'pinp inh')
plot([0 25],[0 25])
xlabel('baseline'); ylabel('laser')
sprintf('%d pinped neurons total',sum(pinped))


figure
plot(baseline(~junk),ev(~junk),'ko')
hold on
plot(baseline(pinped &GT'==1),ev(pinped&GT'==1),'ro');hold on
plot(baseline(pinped &GT'==2),ev(pinped&GT'==2),'bo');hold on
plot(baseline(pinped &GT'==3),ev(pinped&GT'==3),'go');hold on
legend('non','N2B','N2A','wt')
plot([0 50],[0 50]); axis([0 50 0 1000])
xlabel('baseline'); ylabel('laser')


%%% # spikes evoked
nbins=8;
ev_spikes = mean(psth_pinp(pinped,54 + (1:nbins)),2)-baseline(pinped);
ev_spikes = nbins*ev_spikes/1000;
figure
hist(ev_spikes)
xlabel('# evoked spikes - pinped');

nbins=8;
ev_spikes = mean(psth_pinp(pinped_resp_grat,54 + (1:nbins)),2)-baseline(pinped_resp_grat);
ev_spikes = nbins*ev_spikes/1000;
h1 = hist(ev_spikes,-0.05:0.02:0.25)/length(ev_spikes);
ev_spikes = mean(~psth_pinp(~pinped_resp_grat,54 + (1:nbins)),2)-baseline(~pinped_resp_grat);
ev_spikes = nbins*ev_spikes/1000;
h2 = hist(ev_spikes,-0.05:0.02:0.25)/length(ev_spikes);

figure
bar((-0.05:0.02:0.25)+0.01,[h1; h2]')


%%% fraction responsive
frac_resp(1) = sum(resp& lyr==4 & GT'==1 & pinped)/sum( lyr==4 & GT'==1 & pinped);
frac_resp(2) = sum(resp& lyr==4 & GT'==3  & pinped)/sum( lyr==4 & GT'==3  & pinped);
frac_resp(3) = sum(resp& lyr==4 & GT'==2 & pinped)/sum( lyr==4 & GT'==2  & pinped);

figure
bar(frac_resp); ylim([0 2]); ylabel('fraction resp >2'); title('lyr 4')
set(gca,'xticklabel',{'2B pinp','wt pinp','2A pinp'}); ylim([0 1])

%%% fraction responsive
frac_resp(1) = sum(resp& lyr==4 & GT'==1 & ~pinped)/sum( lyr==4 & GT'==1 & ~pinped);
frac_resp(2) = sum(resp& lyr==4 & GT'==3  & ~pinped)/sum( lyr==4 & GT'==3  & ~pinped);
frac_resp(3) = sum(resp& lyr==4 & GT'==2 & ~pinped)/sum( lyr==4 & GT'==2  & ~pinped);

% figure
% bar(frac_resp); ylim([0 2]); ylabel('fraction resp >2'); title('lyr 4')
% set(gca,'xticklabel',{'2B non pinp','wt non pinp','2A pinp'}); ylim([0 1])
%
% %%% fraction responsive all layers
% frac_resp(1) = sum(resp & GT'==1 & ~pinped)/sum(  GT'==1 & ~pinped);
% frac_resp(2) = sum(resp& GT'==3  & ~pinped)/sum(  GT'==3  & ~pinped);
% frac_resp(3) = sum(resp  & GT'==2 & ~pinped)/sum(  GT'==2  & ~pinped);
%
% figure
% bar(frac_resp); ylim([0 2]); ylabel('fraction resp >2'); title('all layer')
% set(gca,'xticklabel',{'2B non pinp','wt non pinp','2A pinp'}); ylim([0 1])

resp=peak(:,2)>=2 & peak(:,1)>=2;

%%wt gain modulation of pinped neurons
peak_run_p=nanmedian(peak(wt &  pinped &  ~junk & resp & ~inh,2))
peak_stat_p=nanmedian(peak(wt  & pinped  & ~junk & resp & ~inh,1))
gain_ind_p=(peak_run_p-peak_stat_p)/(peak_run_p+peak_stat_p)

sprintf('gain wt pinped %f',mean(gain_ind_p))

peak_run=nanmedian(peak(wt   & ~pinped & resp & exc,2))
peak_stat=nanmedian(peak(wt  & ~pinped & resp & exc,1))
gain_ind=(peak_run-peak_stat)/(peak_run+peak_stat)
sprintf('gain wt not pinped %f',mean(gain_ind))

figure
bar([peak_stat_p peak_run_p; peak_stat peak_run])
title('peak resp wt')
legend('stationary','run');
set(gca,'xticklabel',{'pinped','non-pinped'})


%%NR2B gain modualtion of pinped neurons

peak_run_p=nanmedian(peak(N2B & pinped & ~junk & resp & ~inh,2))
peak_stat_p=nanmedian(peak(N2B & pinped & ~junk & resp & ~inh,1))
gain_ind_p=(peak_run_p-peak_stat_p)/(peak_run_p+peak_stat_p)

sprintf('gain n2b pinped %f',mean(gain_ind_p))

peak_run=nanmedian(peak(N2B & ~pinped & resp & exc,2))
peak_stat=nanmedian(peak(N2B & ~pinped & resp & exc,1))
gain_ind=(peak_run-peak_stat)/(peak_run+peak_stat)
sprintf('gain n2b not pinped %f',mean(gain_ind))

figure
bar([peak_stat_p peak_run_p; peak_stat peak_run])
title('peak resp n2b ')
legend('stationary','run');
set(gca,'xticklabel',{'pinped','non-pinped'})

%%NR2A gain modualtion of pinped neurons

peak_run_p=nanmedian(peak(N2A & pinped & resp & exc,2))
peak_stat_p=nanmedian(peak(N2A & pinped & resp & exc,1))
gain_ind_p=(peak_run_p-peak_stat_p)/(peak_run_p+peak_stat_p)

sprintf('gain n2a pinped %f',mean(gain_ind_p))

peak_run=nanmedian(peak(N2A &  ~pinped & resp & exc,2))
peak_stat=nanmedian(peak(N2A &  ~pinped & resp & exc,1))
gain_ind=(peak_run-peak_stat)/(peak_run+peak_stat)
sprintf('gain n2a not pinped %f',mean(gain_ind))

figure
bar([peak_stat_p peak_run_p; peak_stat peak_run])
title('peak resp n2a ')
legend('stationary','run');
set(gca,'xticklabel',{'pinped','non-pinped'})


layer_GT_pinp_bar_ratio(peak(:,1),peak(:,2),GT,lyr,pinped,resp,{'run','stat'});

%%% plot grating response data

plotPinpData(driftspont,wt,lyr<=4,pinped,inh,resp)
plot([0 20],[0 20])
title('spont')

barPinp(driftspont(:,1),wt,N2A,N2B,lyr,pinped,inh,resp)
ylabel('spont')

plotPinpData(peak,wt,lyr<=4,pinped,inh,resp)
plot([0 20],[0 20])
title('peak')

barPinp(peak(:,1),wt,N2A,N2B, lyr,pinped,inh,resp)
ylabel('peak')


plotPinpData(driftF1F0,wt,lyr<=4,pinped,inh,resp)
plot([0 2],[0 2])
title('F1F0')

barPinp(driftF1F0(:,1),wt,N2A,N2B,lyr,pinped,inh,resp)
ylabel('F1F0')


plotPinpData(cvOSI,N2B,lyr<=4,pinped,inh,resp)
plot([0 1],[0 1])
title('cvOSI')

barPinp(cvOSI(:,1),wt,N2A,N2B,lyr,pinped,inh,resp)
ylabel('cvOSI')

plotPinpData(OSI,N2A,lyr<=4,pinped,inh,resp_both)
plot([0 1],[0 1])
title('OSI')

barPinp(OSI(:,1),wt,N2A,N2B,lyr,pinped,inh,resp)
ylabel('OSI')

plotPinpData(cvDSI,wt,lyr==4,pinped,inh,resp)
plot([0 1],[0 1])
title('cv DSI')

barPinp(cvDSI(:,1),wt,N2B,N2A,lyr,pinped,inh,resp)
ylabel('cvDSI')

barPinp(DSI(:,1),wt,N2A,N2B,lyr,pinped,inh,resp)
ylabel('DSI')

lin=STA_exp_var'>0.35 & resp;
barPinp(STA_nx,wt,N2A,N2B,lyr,pinped,inh,resp)
ylabel('nx')

%%cardinal orientation preference ratios
drift_theta_1=size(drift_theta);
drift_theta_1=(drift_theta*180)/pi;
drift_theta_1(drift_theta_1>330)=0;

layerAgePlot_pref_Orient(drift_theta_1(:,1),GT,lyr,inh,pinped,resp ,{'pref Orient' 'Prct total'},'Prefered Orientation Stationary');


%%saptial frequency analysis
driftwpref(driftwpref==0) = 0.005; %tranform SF pref = 0 to 0.005 in order to plot on log scale
driftwpref_1=size(driftwpref);
driftwpref_1=log2(driftwpref);
driftwpref_1=abs(driftwpref_1);


SF_pref_stat =  driftwpref(:,1) >0.006 & driftwpref(:,1)<=0.40 & resp;

barPinp(driftwpref(:,1),wt,N2A,N2B,lyr,pinped,inh,resp)
ylabel('SFpref')

resp = peak(:,1)>=1;

figure
[f,x]=hist(driftwpref(N2B & pinped & resp),0.02:0.04:0.3);
H1=bar(x,f/sum(f));
title 'SF pref N2B'
hold on
% figure
[f1,x1]=hist(driftwpref(wt & pinped & resp),0.02:0.04:0.3);
H2=bar(x1,f1/sum(f1),'b');
ch=get(H2,'child');
set(ch,'facea',.5)
hold on

[f2,x2]=hist(driftwpref(N2A & pinped & resp),0.02:0.04:0.3);
H3=bar(x2,f2/sum(f2),'r');
ch=get(H2,'child');
set(ch,'facea',.5)

f=f/sum(f)
f1=f1/sum(f1)
f2=f2/sum(f2)
h=[f; f1; f2]

figure
bar(x,h')

sf_wt=driftwpref(wt & pinped & resp );
sf_N2B=driftwpref(N2B & pinped & resp );
sf_N2A=driftwpref(N2A & pinped & resp );

kstest2(sf_wt, sf_N2B)
[h p] = kstest2(sf_wt,sf_N2B)

clear f x f1 x1 f2 x2

figure
[f,x]=hist(abs(STA_nx(N2B & pinped & resp)),0.1:0.04:1);
H1=bar(x,f/sum(f),'g');
title 'F1F0 ratio N2B'
hold on


[f1,x1]=hist(abs(STA_nx(wt & pinped & resp)),0.1:0.04:1);
H2=bar(x1,f1/sum(f1),'b');
ch=get(H2,'child');
set(ch,'facea',.5)
hold on

[f2,x2]=hist(abs(STA_nx(N2A & pinped & resp)),0.1:0.04:1);
H2=bar(x1,f1/sum(f1),'r');
ch=get(H2,'child');
set(ch,'facea',.5)
hold on


f=f/sum(f)
f1=f1/sum(f1)
f2=f2/sum(f2)
h=[f; f1; f2]

figure
bar(x,h')
%% simple cell analysis
OS=resp & OSI(:,1)>=0.5;
layerAgePlot_frac_simple(driftF1F0(:,1),GT,lyr,inh, pinped,simple,'F1F0');

figure
[f,x]=hist(driftF1F0(N2B & pinped & OS),0:0.20:2);
H1=bar(x,f/sum(f),'g');
title 'F1F0 ratio N2B'
hold on

% figure
[f1,x1]=hist(driftF1F0(wt & pinped & OS),0:0.20:2);
H2=bar(x1,f1/sum(f1),'b');
ch=get(H2,'child');
set(ch,'facea',.5)
hold on

% [f1,x1]=hist(driftF1F0(wt & ~pinped & OS),0:0.20:2);
% H2=bar(x1,f1/sum(f1),'k');
% ch=get(H2,'child');
% set(ch,'facea',.5)
% hold on

[f2,x2]=hist(driftF1F0(N2A & pinped & OS),0:0.20:2);
H3=bar(x2,f2/sum(f2),'r');
ch=get(H2,'child');
set(ch,'facea',.5)

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
plot(wvform(exc & ~junk,:)','b');
hold on
plot(wvform(inh& ~junk,:)','r'); hold on
plot(wvform(pinped & ~junk,:)','g')

figure
plot(wvform(pinped & ~junk,:)','g')
title('pinped')

l=find(junk)
l=find(N2B & pinped & SF_pref_stat);

%generate STA for pinped cells

for j=1:length(STA_peak)
    
    if ismember(j,l)
        figure
        if ~isempty(STA_peak{1,j})
            colormap jet
            imagesc(STA_peak{1,j}',[-64 64]);
        end
    end
end


driftF1F0(629,1)
OSI(629,1)
DSI(629,1)
driftwpref(629,1)


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

h=session(find(pinped & wt))
% subplot(1,2,2);
% pie([sum(~wt&pinped&~inh) sum(~wt&pinped&inh)],{'broad','narrow'})
% title('2A/2B ko pinped')

%%                       Receptive field STA data

%%%Nx data
clear f x f2 x1
good_STA =STA_exp_var'>=0.1 & STA_ny'<0.39 & responsive_stat;
responsive_stat = peak(:,1)>=2;  % firing rate (responsiveness) criteria for whether cells enter subsequent statistical analysis


barPinp(abs(STA_nx'),wt,N2A,N2B,lyr<=4,pinped,inh, good_STA);
layerAgePlot(abs(STA_ny),age,lyr,inh,good_STA','SF_pref tuning width Stationary ny');
%barPinp(driftspont(:,2),wt,N2A,N2B,lyr<=4,pinped,inh,resp)
ylabel('spont')

STA_nx=STA_nx';
STA_ny=STA_ny';


layerAgeNxNy(abs(STA_nx_j),abs(STA_ny),age,lyr,inh,good_STA,{'nx','ny'},'nx vs ny');




