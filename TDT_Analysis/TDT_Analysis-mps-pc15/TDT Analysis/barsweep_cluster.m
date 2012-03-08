% Matlab codes for reading from TTank
% and plot histgrams and rasters with each orientation
% Jianhua Cang 06-27-03

% connect to the tank and read the block

%close all
clear amp width x0 baseline
[fname, pname] = uigetfile('*.mat','cluster data');
load(fullfile(pname,fname));
% Tank_Name='032806'
% Block_Name='bars1'
TTX = openTTX(Tank_Name,Block_Name); % to initialize



Event_Name_Snip='Snip'
Event_Name_Wave='PDec'
Sample_Interval=0.04096 % 24414.0625Hz
Sample_Number_Snip=64
Dec_Factor=32; % 
plot_duration=6; %in second
hist_int = 0.1;
hist_range=[0:hist_int:9];
axis_range=[0 plot_duration 0 20];
max_events=10000;
deg_per_sec=25;
bar_width = 5;
bar_width_time = bar_width/deg_per_sec;

tetrode_linear=0;

if tetrode_linear
    ch_map = [14 8 10 4 13 7 9 3 11 1 15 5 12 2 16 6];
else
    ch_map = 1:16;
end


%cells = [7 2 ; 9 3; 5 2 ];
cells = [9 2; 10 3; 11 5; 11 6; 12 3 ; 12 5; 13 4; 14 3; 15 1] 
cells = [1 4; 1 6; 5 5; 9 1; 9 2; 9 6; 9 7; 13 2; 13 6]
    
invoke(TTX,'CreateEpocIndexing');

 times = event_times_all(:);
 times = times(times>0);
 time1 = min(times);
 time2= max(times);
MyEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', time1,time2, 1000);
Select_Duration(1) = MyEpocs(2,1); % to exclude events before the first trigger

Screen_Plot_y = [axis_range(3):0.1:axis_range(4)]; % to mark the bar appearing and disappearing time


for cell_n = 1:size(cells,1)
    channel_no = cells(cell_n,1)
    clust_no = cells(cell_n,2)
    hist_fig = figure;
    rast_fig = figure;

    for orientation =0:3
        clear Epocs trial_no Epocs_TS N Spike_TS index TS_xTrg Spike_Timing;
        Epocs=find(MyEpocs(1,:)==orientation+1); 
        trial_no = length(Epocs);
        Epocs_TS = MyEpocs(2,Epocs);

        invoke(TTX,'ResetFilters');
        ecode= invoke(TTX, 'StringToEvCode', 'xTrg'); % to convert string xTrig to code
        invoke(TTX,'SetFilter', ecode, 69, orientation+1, 0); 
                                     % 69 means equal to: 'xTrg=orientation', the last parameter not used for '69'
        N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_Snip, ch_map(channel_no), 0, ...
                    max(time1,Epocs_TS(1)),time2, 'FILTERED');
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
        
        title_text=['Orientation: ' num2str(orientation*90)];

        %raster plot
        figure(rast_fig);
        subplot(2,2,orientation+1); hold on; set(gca, 'yDir','reverse'); 
        axis([0 plot_duration 0 trial_no+1]); title (title_text);
        plot (Spike_Timing, index, '.k', 'MarkerSize',4);

       %% histograms
        figure(hist_fig);

        subplot(2,2,orientation+1);
        bar(hist_range, hist(Spike_Timing, hist_range)/(hist_int*trial_no));  hold on;
        axis(axis_range); title (title_text);
         Bar_Time = episostim_params (orientation*90); % Bar appear and disappear time
         plot(Bar_Time(1)*ones(1, length(Screen_Plot_y)), Screen_Plot_y, 'r-'); % Bar appear
         plot(Bar_Time(2)*ones(1, length(Screen_Plot_y)), Screen_Plot_y, 'r-'); % Bar disappear
        %%%% curve fitting
    
     fit_range = 0:0.1:Bar_Time(2)+0.5;
     Spike_Timing = Spike_Timing(find((Spike_Timing>min(fit_range))&(Spike_Timing<max(fit_range))));

%    obs = hist(Spike_Timing, hist_range)/(hist_int*trial_no);
    %fit_range = hist_range;
    fit_int = fit_range(2)-fit_range(1);
    obs = hist(Spike_Timing, fit_range)/(fit_int*trial_no);
    [min(obs) max(obs) fit_range(find(max(obs))) Bar_Time(2)/5];
    peak_guess = median(fit_range(find(obs> 0.75*max(obs))))
    fit_coeff = nlinfit(fit_range,obs,@rf_fit,[min(obs) max(obs) peak_guess Bar_Time(2)/10])
    baseline(cell_n,orientation+1) = fit_coeff(1);
    amp(cell_n,orientation+1) = fit_coeff(2);
   if orientation<2   
    x0(cell_n,orientation+1) = fit_coeff(3)-Bar_Time(1) + bar_width_time/2;
   else
    x0(cell_n,orientation+1) = Bar_Time(2)-fit_coeff(3) - bar_width_time/2;
   end
    width(cell_n,orientation+1) = abs(fit_coeff(4));
    
    if abs(fit_coeff(4)>(Bar_Time(2)-Bar_Time(1))) | (fit_coeff(2)<0) | ((fit_coeff(4)/fit_coeff(2))>1)
        amp(cell_n,orientation+1)=0;
        x0(cell_n,orientation+1) = 0;
        width(cell_n,orientation+1) = 0;
        baseline(cell_n,orientation+1) = mean(obs);
    end
    
    hold on
    plot(fit_range, rf_fit(fit_coeff,fit_range),'g','LineWidth',1.5);
      
    end

    title_text = sprintf('channel %d cluster %d',channel_no, clust_no);
    figure(rast_fig)
    title(title_text);
    figure(hist_fig);
    title(title_text);
    xlabel('secs');
    ylabel('Hz');
    
    saveas(rast_fig,fullfile(pname,sprintf('rast%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
    saveas(hist_fig,fullfile(pname,sprintf('hist%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
    
end
    baseline
    amp
    x0=x0*deg_per_sec
    width=width*deg_per_sec
    save(fullfile(pname,sprintf('analysis%s_%s',Tank_Name,Block_Name)));
    
invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');