

%function [laserT laserTTL tstamp LFPall] = analyzeEEG_chronux_laser;
% Matlab codes for reading EEG data
% For each condition, goes by epoc and plots power spectrum
% cmn 06-06, based on code by Jianhua Cang 06-27-03



%%% read in cluster data, (only needed to get block name and start/stop time)
%%% then connect to the tank and read the block
clear all
clusterfile=0; %% use cluster file to get info, or not

if clusterfile
    [fname, pname] = uigetfile('*.mat','cluster data');
    load(fullfile(pname,fname));
    times = event_times_all(:);
    times = times(times>0);
    time1 = min(times);
    time2= max(times);
    block = input('Block number : ')
    Block_Name = char(Block_Name(block));
else
    % Tank_Name='06222011_mlr_stim_awake'
    % Block_Name='wn3_72hz'
    pname = uigetdir('C:\data\TDT','block data')
    delims = strfind(pname,'\');
    selected_path = pname(1 :delims(length(delims))-1)
    Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
    Block_Name = pname(delims(length(delims))+1 :length(pname))
    time1= 0;
    time2=0;
    max_events = 50000;
end


%[tsamp vsmooth groomstate] =getBlockVelocity_groom(Tank_Name,Block_Name);

TTX = openTTX(Tank_Name,Block_Name); % to initialize
[laserT laserTTL] = read_laser(TTX);

smoothwindow_secs = 1;
dt = laserT(2)-laserT(1);
smoothwindow = ceil(smoothwindow_secs/dt)

lasersmooth = zeros(size(laserTTL));
for i = smoothwindow+1:length(laserTTL);
    lasersmooth(i)= mean(laserTTL(i-smoothwindow:i));
end
lasersmooth = lasersmooth/5;


Event_Name_EEG='Lfpx'

%%% set time based on first and last timepoints in clustered data
invoke(TTX,'CreateEpocIndexing');

MyEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', time1,time2, 1000);
nEpocs = size(MyEpocs,2);

%%% read in EEG wave info
% N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_EEG, 10, 0, ...
%             MyEpocs(2,1),MyEpocs(2,nEpocs),'All');



LFPall=0;
for tet = 1:4
    
    for ch =1:4
        N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_EEG,4*(tet-1)+ch, 0, ...
            0,0,'All');
        %            N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_EEG,tet, 0, ...
        %             0,0,'All');
        if (N==max_events)
            warning('max number of events aquired');
        end
        W = invoke(TTX, 'ParseEvV', 0, N);
        Wave_TS = invoke(TTX, 'ParseEvInfoV', 0, N, 6);  %   6  = Time Stamp
        W=W(:);
        
        if LFPall==0;
            LFPall=zeros(16, length(W));
        end
        LFPall(4*(tet-1)+ch,:) = W';
        
        
        Wave_TS = Wave_TS(:);
        
        dt=mean(diff(Wave_TS))/32;
        
        tstamp = zeros(length(W),1);
        for i = 1:length(Wave_TS);
            tstamp((i-1)*32+1:i*32) = Wave_TS(i)+(0:31)*dt;
        end
        tstamp = tstamp(1:length(W));
        
        params.Fs = 1/dt;
        params.tapers = [3 5];
        params.fpass = [0 90];
        
        [S t f] = mtspecgramc(W,[3 1],params);
        
        Snorm = zeros(size(S));
        for i = 1:size(S,1);
            Snorm(i,:)=(S(i,:)).*((1:size(S,2)).^1);
        end
        df = max(f)/length(f)
        
        
        %     figure
        %     imagesc(Snorm',[0 10^-7]);
        %     axis xy
        %     set(gca,'YTick',(10:10:80)/df);
        %     set(gca,'YTickLabel',{'10','20','30','40','50','60','70','80'})
        
        
        figure
        %imagesc(Snorm,[0 0.1]);
        imagesc(Snorm',[0 6*10^-8]);
        axis xy
        set(gca,'YTick',(10:10:80)/df);
        set(gca,'YTickLabel',{'10','20','30','40','50','60','70','80'})
        hold on
        plot(laserT,(lasersmooth*100-40),'g');
        axis([0 max(laserT) -40 80/df])
        
        %     figure
        
        %     plot(theta);
        %     hold on
        %     plot(mean(Snorm(:,ceil(12/df):ceil(18/df)),2)+1,'g');
        %     plot(mean(Snorm(:,ceil(20/df):ceil(30/df)),2)+2,'r');
        
        %     plot(gamma+3,'k');
        theta = mean(Snorm(:,ceil(7/df):ceil(10/df)),2);
        gamma = mean(Snorm(:,ceil(55/df):ceil(65/df)),2);
        
        laser_interp = interp1(laserT,lasersmooth,t);
        
        
        %     figure
        %     plot(v_interp,gamma(t),'o');
        %     figure
        %     plot(v_interp,theta(t),'o');
        Smean = mean(Snorm,2)';
        stationary = find(laser_interp<0.02 & Smean<(5*median(Smean)));
        moving = find(laser_interp>0.02 & Smean<(5*median(Smean)));
        
        
        figure
        plot(mean(Snorm(stationary,:),1));
        hold on
        plot(mean(Snorm(moving,:),1),'g');
        
        set(gca,'XTick',(10:10:80)/df);
        set(gca,'XTickLabel',{'10','20','30','40','50','60','70','80'})
        
        title(sprintf('site %d',4*(tet-1)+ch));
        
        if ch==1
            allspectra=figure;
        else
            figure(allspectra);
        end
        subplot(2,2,ch);
        plot(mean(Snorm(stationary,:),1));
        hold on
        plot(mean(Snorm(moving,:),1),'g');
        set(gca,'XTick',(20:20:80)/df);
        set(gca,'XTickLabel',{'20','40','60','80'})
        set(gca,'Xlim',[0 90/df]);
        title(sprintf('channels %d ',4*(tet-1)+ch));
        %         figure
        %         plot(mean(Snorm,2))
        %         title(sprintf('site %d',tet));
    end  %%ch
    
end %% tet

invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');