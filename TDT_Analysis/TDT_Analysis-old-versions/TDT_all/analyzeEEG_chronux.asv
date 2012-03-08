%function analyzeEEG
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
    Tank_Name = '022009_awake_linear';
    Block_Name = 'bars16d1';
    time1= 0;
    time2=0;
    max_events = 50000;
end


[tsamp vsmooth groomstate] =getBlockVelocity_groom(Tank_Name,Block_Name);

TTX = openTTX(Tank_Name,Block_Name); % to initialize

Event_Name_EEG='Lfpx'

%%% set time based on first and last timepoints in clustered data
invoke(TTX,'CreateEpocIndexing');

MyEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', time1,time2, 1000);
nEpocs = size(MyEpocs,2);

%%% read in EEG wave info
% N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_EEG, 10, 0, ...
%             MyEpocs(2,1),MyEpocs(2,nEpocs),'All');
for tet = 13:14
    W=0;
    for ch =1:1
%         N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_EEG,4*(tet-1)+ch, 0, ...
%             0,0,'All');
           N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_EEG,tet, 0, ...
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
    %W= W(:)/4;
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
    imagesc(Snorm',[0 1.5*10^-7]);
    axis xy
    set(gca,'YTick',(10:10:80)/df);
    set(gca,'YTickLabel',{'10','20','30','40','50','60','70','80'})
    hold on
    plot(tsamp,(vsmooth/1.3-40),'g');
    axis([0 max(tsamp) -40 80/df])
    
    %     figure
    
    %     plot(theta);
    %     hold on
    %     plot(mean(Snorm(:,ceil(12/df):ceil(18/df)),2)+1,'g');
    %     plot(mean(Snorm(:,ceil(20/df):ceil(30/df)),2)+2,'r');
    
    %     plot(gamma+3,'k');
    theta = mean(Snorm(:,ceil(7/df):ceil(10/df)),2);
    gamma = mean(Snorm(:,ceil(55/df):ceil(65/df)),2);

    v_interp = interp1(tsamp,vsmooth,t);
    g_interp = interp1(tsamp,groomstate,t);
    %     figure
    %     plot(v_interp,gamma(t),'o');
    %     figure
    %     plot(v_interp,theta(t),'o');
    Smean = mean(Snorm,2)';
    stationary = find(v_interp<1 &g_interp<0.5 & Smean<(5*median(Smean)));
    moving = find(v_interp>5 &g_interp<0.5 & Smean<(5*median(Smean)));
    grooming = find(g_interp>0.5);
    figure
    plot(mean(Snorm(stationary,:),1));
    hold on
    plot(mean(Snorm(moving,:),1),'g');
    if ~isempty(grooming);
        plot(mean(Snorm(grooming,:),1),'r');
    end
    set(gca,'XTick',(10:10:80)/df);
    set(gca,'XTickLabel',{'10','20','30','40','50','60','70','80'})
    
    title(sprintf('site %d',tet));
    figure
    plot(mean(Snorm,2))
     title(sprintf('site %d',tet));
end %% tet


moving30percentbars = mean(Snorm(moving,:),1);

saveas(gcf,fullfile(pname,sprintf('EEG%s',Block_Name)),'fig');

invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');