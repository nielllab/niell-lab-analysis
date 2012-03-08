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
    Tank_Name = '022309_awake_tet';
    Block_Name = 'bars16d1c';
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
for t = 1:4
    W=0;
    for ch =1:4
        N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_EEG,4*(t-1)+ch, 0, ...
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
    onesec=381;
    duration = 4;
    for i = 1:floor(length(W)/onesec)-duration;
        window(i,:) = W((i-1)*onesec+1:(i+duration-1)*onesec);
    end

    blk_size = duration*onesec; %%% approximately 'duration' secs
    nblk = floor(size(W,1)/blk_size)

    W_trim = W(1:blk_size*nblk);
    W_trim = reshape(W_trim,blk_size,nblk);

    W_trim = window';
    %%% take FFT and diplay it
    EEG_power = fft(W_trim,[],1);
    compression_factor = 2
    EEG_condense = condenseData(abs(EEG_power),compression_factor);
    freq = 32/(Wave_TS(5)-Wave_TS(4));
    freq_int = compression_factor*freq/blk_size   %%% this gives frequency scale for x-axis, not sure how to relabel the figure in imagesc
    % figure
    % %imagesc(abs(EEG_condense(1:(180/freq_int),:))',[0 2*10*10^-3]);
    % imagesc(abs(EEG_condense(1:(180/freq_int),:))');

    EEG = abs(EEG_condense(1:(90/freq_int),:))';
    EEG_norm = zeros(size(EEG));
    size(EEG_norm)
    for i = 1:size(EEG,1);
        EEG_norm(i,:) = EEG(i,:).*(1:size(EEG,2));
    end
    figure
    %imagesc(EEG_norm,[0 0.1]);
    imagesc(EEG_norm');
    axis xy
    set(gca,'YTick',(10:10:80)/freq_int);
    set(gca,'YTickLabel',{'10','20','30','40','50','60','70','80'})
hold on
    plot(tsamp,vsmooth/3-15,'g');
    axis([0 max(tsamp) -15 80/freq_int])

    spectrum = mean(EEG_condense(1:(70/freq_int),:),2);
%     figure
%     plot(spectrum.*(1:size(spectrum,1))','r');
%     set(gca,'XTick',(10:10:80)/freq_int);
%     set(gca,'XTickLabel',{'10','20','30','40','50','60','70','80'})

%     figure
    
%     plot(theta);
%     hold on
%     plot(mean(EEG_norm(:,ceil(12/freq_int):ceil(18/freq_int)),2)+1,'g');
%     plot(mean(EEG_norm(:,ceil(20/freq_int):ceil(30/freq_int)),2)+2,'r');

%     plot(gamma+3,'k');
 theta = mean(EEG_norm(:,ceil(7/freq_int):ceil(10/freq_int)),2);
     gamma = mean(EEG_norm(:,ceil(55/freq_int):ceil(65/freq_int)),2); 
t = 1:length(theta);
    v_interp = interp1(tsamp,vsmooth,t-0.5);
  g_interp = interp1(tsamp,groomstate,t-0.5);
%     figure
%     plot(v_interp,gamma(t),'o');
%     figure
%     plot(v_interp,theta(t),'o');
    stationary = find(v_interp<1 );
    moving = find(v_interp>5 );
    grooming = find(g_interp>0.5);
    figure
    plot(mean(EEG_norm(stationary,:),1));
    hold on
    plot(mean(EEG_norm(moving,:),1),'g');
    if ~isempty(grooming);
        plot(mean(EEG_norm(grooming,:),1),'r');
    end
    set(gca,'XTick',(10:10:80)/freq_int);
    set(gca,'XTickLabel',{'10','20','30','40','50','60','70','80'})
    
end %% tet


moving30percentbars = mean(EEG_norm(moving,:),1);

saveas(gcf,fullfile(pname,sprintf('EEG%s',Block_Name)),'fig');

invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');