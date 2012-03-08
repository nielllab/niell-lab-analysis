function barsweep16d_all

% Matlab codes for reading from TTank for sweeping bars in 8 orientations
% plots histgrams and rasters and fits data to a gaussian peak
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03



%%% read in cluster data, then connect to the tank and read the block

    pname = uigetdir('C:\data\','block data')
    delims = strfind(pname,'\');
    selected_path = pname(1 :delims(length(delims))-1)
    Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
    Block_Name = pname(delims(length(delims))+1 :length(pname))

% pname = uigetdir(selected_path,'block data')
% 
% if pname==0;
%     return;
% else
%     nblock=nblock+1;
%     delims = strfind(pname,'\');
%     selected_path = pname(1 :delims(length(delims))-1)
%     Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
%     Block_Name{nblock} = pname(delims(length(delims))+1 :length(pname))
% end

TTX = openTTX(Tank_Name,Block_Name); % to initialize

Event_Name_Snip='Snip'
Sample_Interval=0.04096 % 24414.0625Hz
Sample_Number_Snip=64
Dec_Factor=32; %
plot_duration=3; %in second
hist_int = 0.1;
hist_range=[0:hist_int:9];
axis_range=[0 plot_duration 0 100];
max_events=10000;
deg_per_sec=30;
bar_width = 5;
bar_width_time = bar_width/deg_per_sec;
stim_duration = 3.015;
bar_orients = 0:22.5:337.5
blank_stim=1;
tetrode_linear=0;

if tetrode_linear
    ch_map = [14 8 10 4 13 7 9 3 11 1 15 5 12 2 16 6];
else
    ch_map = 1:16;
end


%   cells = [1 1; 2 1; 3 1; 4 1; 5 1 ; 6 1; 7 1; 8 1; 9 1; 10 1; 11 1; 12 1; 13 1; 14 1; 15 1; 16 1];

cells = [1 1; 5 1; 9 1; 13 1];
%%% set time based on first and last timepoints in clustered data
invoke(TTX,'CreateEpocIndexing');




% if exist('blockID')
%     times = event_times_all(find(blockID==block));
% else
%     times = event_times_all(:);
% end
%  times = times(times>0);
trgname = 'xTrg';

MyEpocs = invoke(TTX, 'GetEpocsV', trgname, 0,0, 1000)


Screen_Plot_y = [axis_range(3):0.1:axis_range(4)]; % to mark the bar appearing and disappearing time


for cell_n = 1:size(cells,1)
    channel_no = cells(cell_n,1)
    clust_no = cells(cell_n,2)
    hist_fig = figure;
    rast_fig = figure;
    
    
    
    for orientation =0:size(bar_orients,2)+blank_stim-1;
        clear Epocs trial_no Epocs_TS N Spike_TS index TS_xTrg Spike_Timing;
        Epocs=find(MyEpocs(1,:)==orientation+1);
        trial_no = length(Epocs);
        Epocs_TS = MyEpocs(2,Epocs);
        
        invoke(TTX,'ResetFilters');
        ecode= invoke(TTX, 'StringToEvCode', trgname); % to convert string xTrig to code
        invoke(TTX,'SetFilter', ecode, 69, orientation+1, 0);
        % 69 means equal to: 'xTrg=orientation', the last parameter not used for '69'
        N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_Snip, ch_map(channel_no), 0, ...
            0,0, 'FILTERED');
        if (N==max_events)
            warning('max number of events aquired');
        end
        Spike_TS = invoke(TTX, 'ParseEvInfoV', 0, N, 6);  %   6  = Time Stamp
        
        
        
        % Now for all events, we find out when the xTrig was:
        index = zeros(1,N);
        TS_xTrg=index;
        for i = 1:size(Epocs_TS,2)-1;
            epochSpikes = find(Spike_TS>=Epocs_TS(i) & Spike_TS<Epocs_TS(i+1));
            index(epochSpikes)=i;
            TS_xTrg(epochSpikes)=Epocs_TS(i);
        end
        Spike_Timing=Spike_TS-TS_xTrg;
        
        
        %%% calculate total spikes
        R(cell_n,orientation+1) = sum(Spike_Timing<stim_duration)/stim_duration;
        
        title_text=['Orientation: ' num2str(orientation*22.5)];
        
        if orientation<16
            %%%raster plot
            position = [6 9 8 7 4 1 2 3 5];
            figure(rast_fig);
            subplot(4,4,orientation+1); hold on; set(gca, 'yDir','reverse');
            axis([0 plot_duration 0 trial_no+1]); title (title_text);
            plot (Spike_Timing, index, '.k', 'MarkerSize',4);
            
            
            %% histograms
            figure(hist_fig);
            subplot(4,4,orientation+1);
            bar(hist_range, hist(Spike_Timing, hist_range)/(hist_int*trial_no));  hold on;
            axis(axis_range); title (title_text);
            Bar_Time = episostim_params (orientation*90); % Bar appear and disappear time
            plot(Bar_Time(1)*ones(1, length(Screen_Plot_y)), Screen_Plot_y, 'r-'); % Bar appear
            plot(Bar_Time(2)*ones(1, length(Screen_Plot_y)), Screen_Plot_y, 'r-'); % Bar disappear
            
            %%%% curve fitting
            
            fit_range = 0:0.1:Bar_Time(2)+0.5;
            Spike_Timing = Spike_Timing(find((Spike_Timing>min(fit_range))&(Spike_Timing<max(fit_range))));
            fit_int = fit_range(2)-fit_range(1);
            obs = hist(Spike_Timing, fit_range)/(fit_int*trial_no);
            [min(obs) max(obs) fit_range(find(max(obs))) Bar_Time(2)/5];
            peak_guess = median(fit_range(find(obs> 0.75*max(obs))))
            fit_coeff = nlinfit(fit_range,obs,@rf_fit,[min(obs) max(obs) peak_guess Bar_Time(2)/10])
            baseline(cell_n,orientation+1) = fit_coeff(1);
            amp(cell_n,orientation+1) = fit_coeff(2);
            if orientation<4
                x0(cell_n,orientation+1) = fit_coeff(3) + bar_width_time/2;
            else
                x0(cell_n,orientation+1) = stim_duration-fit_coeff(3) - bar_width_time/2;
            end
            width(cell_n,orientation+1) = abs(fit_coeff(4));
            
            %%% look for aberrant results, and set all values to zero
            if abs(fit_coeff(4)>(Bar_Time(2)-Bar_Time(1))) | (fit_coeff(2)<0) | ((fit_coeff(4)/fit_coeff(2))>1)
                amp(cell_n,orientation+1)=0;
                x0(cell_n,orientation+1) = 0;
                width(cell_n,orientation+1) = 0;
                baseline(cell_n,orientation+1) = mean(obs);
            end
            
            
            hold on
            plot(fit_range, rf_fit(fit_coeff,fit_range),'g','LineWidth',1.5);
            
        end
    end   %orientation
    
    title_text = sprintf('channel %d cluster %d',channel_no, clust_no);
    figure(rast_fig)
    title(title_text);
    figure(hist_fig);
    title(title_text);
    xlabel('secs');
    ylabel('Hz');
    
    %     saveas(rast_fig,fullfile(pname,sprintf('bar_rast%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
    %     saveas(hist_fig,fullfile(pname,sprintf('bar_hist%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
    %
end  %%% cell
baseline
amp
x0=x0*deg_per_sec
width=width*deg_per_sec

% if use_afile
%     barsweep_R = R;
%     bars_x0 = x0;
%     bars_peakwidth = width;
%     bars_amp = amp;
%     bars_baseline = baseline;
%     save(afile, 'barsweep_R', 'bars_x0', 'bars_peakwidth', 'bars_amp', 'bars_baseline', 'bar_orients','-append');
% end
% clear event_times_all etimes_old idx_all tm used score c_score channel_times
% save(fullfile(pname,sprintf('analysis%s_%s',Tank_Name,Block_Name)));

invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');