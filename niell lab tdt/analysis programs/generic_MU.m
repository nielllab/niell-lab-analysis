function generic_MU

% Matlab codes for reading from TTank for sweeping bars in 8 orientations
% plots histgrams and rasters and fits data to a gaussian peak
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03

clear all

%%% read in cluster data, then connect to the tank and read the block
pname = uigetdir('C:\data\','block data')
delims = strfind(pname,'\');
selected_path = pname(1 :delims(length(delims))-1)
Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
Block_Name = pname(delims(length(delims))+1 :length(pname))

stim_duration = input('duration : ');
nrows = input('rows : ');
ncols = input('cols : ');
panels = input('panels : ');


plot_duration=stim_duration; %in second
hist_int = 0.1;
hist_range=[0:hist_int:9];
axis_range=[0 plot_duration 0 50];
deg_per_sec=30;

bar_orients = 0:22.5:337.5
blank_stim=1;

chans = [1 5 9 13];

flags = struct('visStim',1,'MUspike',1);
data = getTDTdata(Tank_Name,Block_Name,chans,flags)
allEpocs = data.stimEpocs;

for ch=chans;
    
    hist_fig = figure;
    for rep = 1:panels
        
        
        rast_fig = figure;
        Spike_TS =data.MUspikeT{ch};
        N= length(Spike_TS);
        
        for c =1:nrows*ncols;
            clear  index TS_xTrg Spike_Timing;
            
            orientation = (c-1)*panels+rep;
            %%% find epocs for this orientation
            condEpocs=find(allEpocs(1,:)==orientation);
            trial_no = length(condEpocs);
            Epocs_TS = allEpocs(2,condEpocs);
            
            
            % Now for all events, we find out when the xTrig was:
            index = zeros(1,N);
            TS_xTrg=index;
            for i = 1:size(Epocs_TS,2)-1;
                epochSpikes = find(Spike_TS>=Epocs_TS(i) & Spike_TS<Epocs_TS(i+1));
                index(epochSpikes)=i;
                TS_xTrg(epochSpikes)=Epocs_TS(i);
            end
            
            %%% subtract off trial start time, and only keep spikes for this
            %%% orientation i.e. index>0
            Spike_Timing=Spike_TS(index>0)-TS_xTrg(index>0);
            index = index(index>0);
            
            title_text=['theta = ' num2str(orientation*22.5)];
       
                %%%raster plot
                position = [6 9 8 7 4 1 2 3 5];
                figure(rast_fig);
                subplot(nrows,ncols,c); hold on; set(gca, 'yDir','reverse');
                axis([0 plot_duration 0 trial_no+1]); 
                plot (Spike_Timing, index, '.k', 'MarkerSize',4);
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
            plot(hist_range, hist(Spike_Timing, hist_range)/(hist_int*trial_no),color);
            hold on;
            axis(axis_range);
             set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])    
        end   %orientation
        
%         title_text = sprintf('channel %d',ch);
%         figure(rast_fig)
%         title(title_text);
%         figure(hist_fig);
%         title(title_text);
%         xlabel('secs');
%         ylabel('Hz');
        
        %     saveas(rast_fig,fullfile(pname,sprintf('bar_rast%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
        %     saveas(hist_fig,fullfile(pname,sprintf('bar_hist%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
        
    end
    %
end  %%% cell


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
