%function analyzeEEG
% Matlab codes for reading EEG data
% For each condition, goes by epoc and plots power spectrum
% cmn 06-06, based on code by Jianhua Cang 06-27-03



%%% read in cluster data, (only needed to get block name and start/stop time)
%%% then connect to the tank and read the block

%%% read in cluster data, then connect to the tank and read the block
clear all
close all
cells =1;
[fname, pname] = uigetfile('*.mat','cluster data');
load(fullfile(pname,fname));
for i =1:length(Block_Name);
    sprintf('%d : %s ',i,Block_Name{i})
end
block = input('which block to analyze ? ');
Block_Name = Block_Name{block}

TTX = openTTX(Tank_Name,Block_Name); % to initialize
time1= 0;
time2=0;
max_events = 50000;

[tsamp vsmooth] = getBlockVelocity(Tank_Name,Block_Name);

thresh_velocity=1;

[afname, apname] = uigetfile('*.mat','analysis data');
if afname~=0
    afile_bar = fullfile(apname,afname)
    save(afile_bar, 'afile_bar','apname','-append'); %%% to make sure it's not overwritten
    load(afile_bar);
    use_afile=1;
else
    use_afile=0;
    %%% select which units to analyze (channel, cluster number)

end


[tsamp vsmooth groomstate] =getBlockVelocity_groom(Tank_Name,Block_Name);

TTX = openTTX(Tank_Name,Block_Name); % to initialize

Event_Name_EEG='Lfpx'


invoke(TTX,'CreateEpocIndexing');

MyEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', time1,time2, 1000);
nEpocs = size(MyEpocs,2);
cell_n=0;
%%% read in EEG wave info
% N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_EEG, 10, 0, ...
%             MyEpocs(2,1),MyEpocs(2,nEpocs),'All');
for tet = 1:4
    %tet=1;
  
   W=0;
    for ch =1:4
        N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_EEG,4*(tet-1)+ch, 0, ...
            0,0,'All');

        if (N==max_events)
            warning('max number of events aquired');
        end
        W_onechan = invoke(TTX, 'ParseEvV', 0, N);
        Wave_TS = invoke(TTX, 'ParseEvInfoV', 0, N, 6);  %   6  = Time Stamp
        W = W+W_onechan(:);
    end%% ch

    %%% reshape wave data so that one FFT block runs down a column
    %%% samples at 381 Hz;

    W= W(:)/4;
    Wave_TS = Wave_TS(:);

    dt=mean(diff(Wave_TS))/32;

    tstamp = zeros(length(W),1);
    for i = 1:length(Wave_TS);
        tstamp((i-1)*32+1:i*32) = Wave_TS(i)+(0:31)*dt;
    end
    tstamp = tstamp(1:length(W));

    [ph amp t] = phaseFromTimeseries(W,tstamp,9,1);

    block_times = event_times_all-(block-1)*10^5;

    tetclusters = find(cells(:,1)==4*(tet-1)+1);
    for i = 1:length(tetclusters);
 cell_n=cell_n+1;
        c=tetclusters(i);
        ch = cells(c,1);
        clust = cells(c,2);

        channel_times =squeeze(block_times(ch,:));
        used = find(idx_all(ch,:)==clust & channel_times>0 & channel_times<10^5);

        spiketimes = channel_times(used);
        spiketimes = spiketimes(spiketimes>min(t) & spiketimes<max(t));
        
        allspikes{cell_n} = spiketimes;
        all_lfp{cell_n} = W;
        all_t_lfp{cell_n} = tstamp;
        all_t_movement{cell_n} = tsamp;
        all_v{cell_n} = vsmooth;
        
        
        v_interp = interp1(tsamp,vsmooth,spiketimes);
        g_interp = interp1(tsamp,groomstate,spiketimes);

          spiket = spiketimes(v_interp>5);
          
        phases = interp1(t, ph, spiket);
        figure
        hist(mod(phases,2*pi));
        title(sprintf('tet %d cluster %d',ch,clust));
        
        resamp_t = 0:.002:max(spiketimes);
        Wsamp = interp1(tstamp,W,resamp_t);
        
        
        spikehist = hist(spiket,resamp_t);
        figure
        plot(abs(fftshift(fft(spikehist))));
               title(sprintf('tet %d cluster %d',ch,clust));
%         figure
%         plot(abs(fftshift(fft(spikehist).*fft(Wsamp))));
%          title(sprintf('tet %d cluster %d',ch,clust));    
         
     
         range=50;
         dt = .005
         for i = -range:range;
             sta(i+range+1)=mean(interp1(tstamp,W,spiket+i*dt));
         end
         figure
         plot((-range:range)*dt*1000,sta);
          title(sprintf('tet %d cluster %d',ch,clust));
         figure
         plot((-range:range) / (range*dt*2),abs(fftshift(fft(sta))))
%           title(sprintf('tet %d cluster %d',ch,clust));
% figure
%         plot(abs(conv(ones(50,1), fft(spikehist).*abs(fft(Wsamp))./fft(Wsamp))));
        
%          figure
%             plot(ifft(fft(spikehist).*fft(Wsamp)));
         
        
%         figure
%         hist(interp1(tstamp,W,spiketimes(v_interp>5)));
%      title(sprintf('tet %d cluster %d',ch,clust));
     
     mean(interp1(tstamp,W,spiketimes(v_interp>5)))
     
     dt=mean(diff(Wave_TS))/32;
     params.tapers = [6 11];
     params.fpass = [0 90];
     params.Fs = 1/dt;
     [C phi S12 S1 S2 win_t f] = cohgramcpt(W,spiketimes',[5 1],params);
     
      df = max(f)/length(f)
     figure
     imagesc(C');
     axis xy
     set(gca,'YTick',(10:10:80)/df);
    set(gca,'YTickLabel',{'10','20','30','40','50','60','70','80'})
     
     v_interp = interp1(tsamp,vsmooth,win_t);
     
     moving = v_interp>2 & ~isnan(C(:,1))';
     
%          figure
%     plot(mean(C(moving,:),1));
%          set(gca,'XTick',(10:10:80)/df);
%     set(gca,'XTickLabel',{'10','20','30','40','50','60','70','80'})
%      
    figure
    plot(abs(mean(C(moving,:).*exp(sqrt(-1)*phi(moving,:)),1)));
         set(gca,'XTick',(10:10:80)/df);
    set(gca,'XTickLabel',{'10','20','30','40','50','60','70','80'})
    end



end %% tet


invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');