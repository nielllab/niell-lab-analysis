function raster_generic(Tank_Name,Block_Name, channel_list, nrows, ncols, bin_size,max_rate)
% Matlab codes for reading from TTank for sweeping bars in 8 orientations
% plots histgrams and rasters and fits data to a gaussian peak
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03


TTX = openTTX(Tank_Name,Block_Name); % to initialize

Event_Name_Snip='Snip'

invoke(TTX,'CreateEpocIndexing');   %%% Epocs are timeranges for conditions

MyEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', 0,0, 1000);

duration = MyEpocs(2,2)-MyEpocs(2,1)  %%% calculate duration of trials by length of first Epoc
max_events = 10000;
for channel_no=channel_list;

    hist_fig = figure;
    rast_fig = figure;

    nConds = max(MyEpocs(1,:))


    for cond =1:min(nConds,nrows*ncols);
        clear Epocs trial_no Epocs_TS N Spike_TS index TS_xTrg Spike_Timing;
        Epocs=find(MyEpocs(1,:)==cond);
        nTrials = length(Epocs);
        Epocs_TS = MyEpocs(2,Epocs);

        invoke(TTX,'ResetFilters');
        ecode= invoke(TTX, 'StringToEvCode', 'xTrg'); % to convert string xTrig to code
        invoke(TTX,'SetFilter', ecode, 69, cond, 0);
        % 69 means equal to: 'xTrg=cond', the last parameter not used for '69'
        N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_Snip, channel_no, 0, ...
            0,0, 'FILTERED');
        if (N==max_events)
            warning('max number of events acquired');
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

        %%%raster plot
    
        figure(rast_fig);
        subplot(nrows,ncols,cond); hold on; set(gca, 'yDir','reverse');
        axis([0 duration 0 nTrials+1]); 
        plot (Spike_Timing, index, '.k', 'MarkerSize',4);

        %% histograms
        figure(hist_fig);
        subplot(nrows, ncols,cond);
        hist_range = 0:bin_size:duration;
        bar(hist_range, hist(Spike_Timing, hist_range)/(bin_size*nTrials));  hold on;
        axis([0 duration 0 max_rate])
    end  %%% cond

    title_text = sprintf('channel %d',channel_no);
    figure(rast_fig)
    title(title_text);
    figure(hist_fig);
    title(title_text);

end  %%% channel_no


invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');