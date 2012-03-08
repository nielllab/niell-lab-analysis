%function movie_cluster
% Matlab codes for reading from TTank
% and plot histgrams and rasters with each orientation
% Jianhua Cang 06-27-03

% connect to the tank and read the block

%close all
clear amp width x0 baseline
[fname, pname] = uigetfile('*.mat','cluster data');
oldpname = pname;
load(fullfile(pname,fname));
pname = oldpname;
TTX = openTTX(Tank_Name,Block_Name); % to initialize

Event_Name_Snip='Snip'
Event_Name_Wave='PDec'
Sample_Interval=0.04096 % 24414.0625Hz
Sample_Number_Snip=64
Dec_Factor=32; % 
plot_duration=121; %in second
hist_int = 0.1;
hist_range=[0:hist_int:plot_duration];
axis_range=[0 plot_duration 0 20];
max_events=50000;


tetrode_linear=0;

if tetrode_linear
    ch_map = [14 8 10 4 13 7 9 3 11 1 15 5 12 2 16 6];
else
    ch_map = 1:16;
end


%cells = [7 2 ; 9 3; 5 2 ];
cells = [2 4; 3 2; 5 3; 7 4 ; 7 7; 8 4; 9 5; 13 4; 14 4; 15 3];
cells = [15 2];

cells = [12 5; 12 3; 15 2];
invoke(TTX,'CreateEpocIndexing');

MyEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', 0,max_time, 1000);
%Select_Duration(1) = MyEpocs(2,1); % to exclude events before the first trigger

Screen_Plot_y = [axis_range(3):0.1:axis_range(4)]; % to mark the bar appearing and disappearing time
  times = event_times_all(:);
 times = times(times>0);
 time1 = min(times);
 time2= max(times);

for cell_n = 1:size(cells,1)
    channel_no = cells(cell_n,1)
    clust_no = cells(cell_n,2)
    hist_fig = figure;
    rast_fig = figure;

    for orientation =0
        clear Epocs trial_no Epocs_TS N Spike_TS index TS_xTrg Spike_Timing;
        Epocs=find(MyEpocs(1,:)==orientation+1); 
        trial_no = length(Epocs);
        Epocs_TS = MyEpocs(2,Epocs);

        invoke(TTX,'ResetFilters');
        ecode= invoke(TTX, 'StringToEvCode', 'xTrg'); % to convert string xTrig to code
        invoke(TTX,'SetFilter', ecode, 69, orientation+1, 0); 
                                     % 69 means equal to: 'xTrg=orientation', the last parameter not used for '69'
        N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_Snip, ch_map(channel_no), 0, ...
                    time1, time2, 'FILTERED');
        if (N==max_events)
            warning('max number of events aquired');
        end
        Spike_TS = invoke(TTX, 'ParseEvInfoV', 0, N, 6);  %   6  = Time Stamp

       channel_times =squeeze(event_times_all(channel_no,:));
        spike_no=zeros(N,1);        
        n=1;
        for i = 1:N        
            while channel_times(n)~=Spike_TS(i)
                n = n+1;
            end
            spike_no(i) = n;
        end        
        clust_idx = idx_all(channel_no,spike_no);


        % Now for all events, we find out when the xTrig was:
        index = zeros(1,N);
        TS_xTrg=index;
        for i = 1:size(Epocs_TS,2)-1;
            epochSpikes = find(Spike_TS>=Epocs_TS(i) & Spike_TS<Epocs_TS(i+1));
            index(epochSpikes)=i;
            TS_xTrg(epochSpikes)=Epocs_TS(i);
        end

        Spike_Timing=Spike_TS-TS_xTrg;

        Spike_Timing = Spike_Timing(find(clust_idx==clust_no));
        index=index(find(clust_idx==clust_no));
        
        title_text=Block_Name;

        %raster plot
        figure(rast_fig);
        hold on; set(gca, 'yDir','reverse'); 
        axis([0 plot_duration 0 trial_no+1]); title (title_text);
        plot (Spike_Timing, index, '.k', 'MarkerSize',2);
        
        rate = zeros(size(hist_range,2),max(index));
        for i = 1:max(index);
            rate(:,i) = hist(Spike_Timing(find(index==i)),hist_range)/(hist_int*trial_no);
            fano = var(rate')./mean(rate');
        end
          cc = corrcoef(rate);
          figure
          imagesc(cc,[0 1]);
        colorbar
        
 
        
       %% histograms
        figure(hist_fig);

       
        rate = hist(Spike_Timing, hist_range)/(hist_int*trial_no);
        bar(hist_range, rate);  hold on;
        axis(axis_range); title (title_text);
       
        figure
        [y x] =hist(rate,0:0.2:10)
        
       bar(x,y)
    end
    title_text = sprintf('channel %d cluster %d',channel_no, clust_no);
    figure(rast_fig)
    title(title_text);
    figure(hist_fig);
    title(title_text);
    xlabel('secs');
    ylabel('Hz');
end

    
invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');