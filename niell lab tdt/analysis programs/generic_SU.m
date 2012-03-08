function generic_SU
% Matlab codes for reading from TTank for sweeping bars in 8 orientations
% plots histgrams and rasters and fits data to a gaussian peak
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03



%%% read in cluster data, then connect to the tank and read the block
clear all
cells =1;
[fname, pname] = uigetfile('*.mat','cluster data');
load(fullfile(pname,fname));
for i =1:length(Block_Name);
    sprintf('%d : %s ',i,Block_Name{i})
end
block = input('which block to analyze ? ');
Block_Name = Block_Name{block}

stim_duration = input('duration : ');
nrows = input('rows : ');
ncols = input('cols : ');
panels = input('panels : ');



deg_per_sec=30;


plot_duration=stim_duration; %in second

hist_int = plot_duration/20;
hist_range=[0:hist_int:plot_duration];
axis_range=[0 plot_duration 0 30];


blank_stim=1;

[afname, apname] = uigetfile('*.mat','analysis data');
if afname~=0
    afile = fullfile(apname,afname);
    load(afile);
    use_afile=1;
else
    use_afile=0;
    %%% select which units to analyze (channel, cluster number)
    
end
cells


%printfig = input('print figures? ');
printfig=0;
for cell_n = 1:size(cells,1)
    % for cell_n=9:9
    cell_n
    channel_no = cells(cell_n,1)
    clust_no = cells(cell_n,2)
    channel_times =spikeT{cell_n} - (block-1)*10^5;
    times = channel_times(channel_times>0 & channel_times<10^5);
    hist_fig = figure('Name',sprintf('unit %d %d',channel_no,clust_no)) 
    for rep =1:panels
        
        rast_fig = figure('Name',sprintf('unit %d %d rep %d',channel_no,clust_no,rep));
        for c =1:nrows*ncols;
            orientation = (c-1)*panels+rep;
            [Spike_Timing index numtrials] = getTrials(stimEpocs{block},times, orientation, stim_duration);
            
            %%% calculate total spikes
            R(cell_n,orientation+1) = sum(Spike_Timing<stim_duration)/(stim_duration*numtrials);
            
            title_text=['Orientation: ' num2str(orientation*22.5)];
            
            
            figure(rast_fig);
            subplot(nrows,ncols,c)
            hold on; set(gca, 'yDir','reverse');
            axis([0 plot_duration 0 numtrials+1]);
            %title (title_text);
            %
            %         plot (Spike_Timing, index, '.k', 'MarkerSize',4);
            
            plot ([Spike_Timing; Spike_Timing], [index-0.25;index+0.25], 'k', 'MarkerSize',4);
            set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])
            
            %% histograms
            figure(hist_fig);
            subplot(nrows,ncols,c);
            hold on
            if rep ==1
                color = 'b';
            else
                color = 'r';
            end
            plot(hist_range, hist(Spike_Timing, hist_range)/(hist_int*numtrials),color);
            hold on;
            axis(axis_range);
             set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])         
%             if orientation==1
%                 title_text = sprintf('ch%d c%d rep%d',channel_no,clust_no, rep);
%             end
            
            
        end %orientation
        saveas(rast_fig,fullfile(pname,sprintf('generic_rast_move%d%s_%d_%d',rep,Block_Name,channel_no,clust_no)),'fig');
        saveas(hist_fig,fullfile(pname,sprintf('generic_hist_move%d%s_%d_%d',rep,Block_Name,channel_no,clust_no)),'fig');
        
        % if use_afile
        %     barsweep_R = R;
        %     bars_x0 = x0;
        %     bars_peakwidth = width;
        %     bars_amp = amp;
        %     bars_baseline = baseline
        %     save(afile, 'barsweep_R', 'bars_x0', 'bars_peakwidth', 'bars_amp', 'bars_baseline', ...
        %         'bar_orients',  'bars_theta', 'bars_OSI', 'bars_A1', 'bars_A2', 'bars_w', 'bars_B','bars_null', 'bars_spont','rf_width','bars_dsi','-append');
        % end
        
    end
end

