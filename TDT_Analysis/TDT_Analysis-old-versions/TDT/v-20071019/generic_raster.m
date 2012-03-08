function generic_raster(Tank_Name,Block_Name,channel_list,nrows,ncols,bin_size,max_rate,varargin)
% Matlab codes for reading from TTank for sweeping bars in 8 orientations
% plots histgrams and rasters and fits data to a gaussian peak
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03

%%% read in cluster data, then connect to the tank and read the block

figurescoexistwithtdt = 0;

isorientationbatch = 0;
if (~isempty(varargin))
    isorientationbatch = varargin{1};
end

% Tank_Name='072607_wt_linear'
% Block_Name='bars16d1b'
% nrows = 4;
% ncols = 4;
% bin_size=0.1;
% max_rate = 20;
% channel_list = 1:16;

TTX = openTTX(Tank_Name,Block_Name); % to initialize

Event_Name_Snip = 'Snip';

%   cells = [1 1; 2 1; 3 1; 4 1; 5 1 ; 6 1; 7 1; 8 1; 9 1; 10 1; 11 1; 12 1; 13 1; 14 1; 15 1; 16 1];

cells = [1 1; 5 1; 9 1; 13 1];

%%% set time based on first and last timepoints in clustered data
invoke(TTX,'CreateEpocIndexing');   %%% Epocs are timeranges for conditions

MyEpocs = invoke(TTX,'GetEpocsV','xTrg',0,0,1000);
if (isnan(MyEpocs))
    error('Cannot find triggers in this block');
end
    
duration = MyEpocs(2,2) - MyEpocs(2,1);  %%% calculate duration of trials by length of first Epoc
max_events = 10000;
numchannels = length(channel_list);
for ch = 1:numchannels
    channel_no = channel_list(ch);
    
    if (numchannels == 1)
        hist_fig = figure(1);
        set(hist_fig,'Units','normalized');
        if figurescoexistwithtdt
            set(hist_fig,'Position',[0.7613    0.3857    0.2381    0.2857]);
            % set(hist_fig,'Position',[1280 406 400 300]);
        else
            set(hist_fig,'Position',[0    0.3857    1    0.2857]);
        end

        rast_fig = figure(2);
        set(rast_fig,'Units','normalized');
        if figurescoexistwithtdt
            set(rast_fig,'Position',[0.7613    0.0295    0.2381    0.2857]);
            % set(rast_fig,'Position',[1280 32 400 300]);
        else
            set(rast_fig,'Position',[0    0.0295    1    0.2857]);
        end
            
    else
        side = ceil(sqrt(numchannels));
        row = floor((ch-1)/side);
        col = mod((ch-1),side);
        hist_fig = subfig((ch-1)*2+1,side*2,side,row*side*2 + col + 1);
        rast_fig = subfig(ch*2,side*2,side,row*side*2 + side + col + 1);
    end

    nConds = max(setdiff(MyEpocs(1,:),[999])); % max condition which is not stop code

    for cond = 1:min(nConds,nrows*ncols);
        clear Epocs trial_no Epocs_TS N Spike_TS index TS_xTrg Spike_Timing;
        Epocs=find(MyEpocs(1,:)==cond);
        nTrials = length(Epocs);
        Epocs_TS = MyEpocs(2,Epocs);

        invoke(TTX,'ResetFilters');
        ecode= invoke(TTX, 'StringToEvCode', 'xTrg'); % to convert string xTrig to code
        invoke(TTX,'SetFilter', ecode, 69, cond, 0);
        % 69 means equal to: 'xTrg=cond', the last parameter not used for '69'
        N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_Snip, channel_no, 0, 0,0,'FILTERED');
        if (N==max_events)
            disp('max number of events acquired');
        end
        Spike_TS = invoke(TTX, 'ParseEvInfoV', 0, N, 6);  %   6  = Time Stamp

        % Now for all events, we find out when the xTrig was:
        index = zeros(1,N);
        TS_xTrg=index;
        for i = 1:length(Epocs_TS);
            epochSpikes = find(Spike_TS>Epocs_TS(i) & Spike_TS<(Epocs_TS(i)+duration));
            index(epochSpikes)=i;
            TS_xTrg(epochSpikes)=Epocs_TS(i);
        end

        Spike_Timing=Spike_TS-TS_xTrg;

        %          %%% calculate total spikes
        %          R(cell_n,orientation+1) = sum(Spike_Timing<stim_duration)/stim_duration;

        if ~isorientationbatch
            condvalue = cond-1;
        else
            orientations = linspace(0,360,nConds+1);
            % orientations = orientations(1:end-1);
            condvalue = round(orientations(cond));
        end

        title_text = sprintf('ch%d cond %d',channel_no,condvalue);

        %% raster plot
        figure(rast_fig);
        subplot(nrows,ncols,cond); hold on; set(gca,'yDir','reverse');
        axis([0 duration 0 nTrials+1]); 
        plot(Spike_Timing, index, '.k', 'MarkerSize', 4);
        title(title_text);

        %% histograms
        figure(hist_fig);
        subplot(nrows, ncols,cond);
        hist_range = 0:bin_size:duration;
        if (nTrials > 0)
            bar(hist_range, hist(Spike_Timing, hist_range) / (bin_size*nTrials));  
        end
        hold on;
        axis([0 duration 0 max_rate])
        title(title_text);
    end  %%% cond

end  %%% channel_no

invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');
