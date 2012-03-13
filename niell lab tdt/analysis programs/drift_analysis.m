function drift_analysis
% Matlab codes for reading from TTank for sweeping bars in 8 orientations
% plots histgrams and rasters and fits data to a gaussian peak
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03



%%% read in cluster data, then connect to the tank and read the block
clear all
SU = menu('recording type','multi-unit','single unit')-1

cells =1;
if SU
    [fname, pname] = uigetfile('*.mat','cluster data');
    load(fullfile(pname,fname));
    block = listdlg('ListString',Block_Name,'SelectionMode','single');
    Block_Name = Block_Name{block}
    [afname, apname] = uigetfile('*.mat','analysis data');
    noisepname = apname;
    afile = fullfile(apname,afname);
    load(afile);
    use_afile=1;
    cells
else
    pname = uigetdir('C:\data\','block data')
    delims = strfind(pname,'\');
    selected_path = pname(1 :delims(length(delims))-1)
    Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
    Block_Name = pname(delims(length(delims))+1 :length(pname))
    nchan = input('number of channels : ');
    flags = struct('visStim',1,'MUspike',1);
    data = getTDTdata(Tank_Name,Block_Name,1:4:nchan,flags);
    
end

stim_duration = input('duration : ');
nrows = input('# orients : ');
ncols = input('# sfs : ');
tf = input('temp freqs : ');
latency = input('latency (0.05]');
panels= length(tf);

plot_duration=stim_duration; %in second

hist_int = plot_duration/20;
hist_range=[0:hist_int:plot_duration];
axis_range=[0 plot_duration 0 30];

if SU
    cell_range = 1:size(cells,1)
else
    cell_range=1:4:nchan;
end

for cell_n = cell_range;
    % for cell_n=9:9
    cell_n
    if SU
        channel_no = cells(cell_n,1)
        clust_no = cells(cell_n,2)
        channel_times =spikeT{cell_n} - (block-1)*10^5;
        times = channel_times(channel_times>0 & channel_times<10^5);
        hist_fig = figure('Name',sprintf('unit %d %d',channel_no,clust_no))
    else
        hist_fig = figure('Name',sprintf('channel %d',cell_n))
    end
    
    %%% spont and flicker
    spontfig=figure
    for rep=1:2
        cond = panels*nrows*ncols+rep;
        
        if SU
            [Spike_Timing index numtrials] = getTrialsSU(stimEpocs{block},times, cond, stim_duration);
        else
            [Spike_Timing index numtrials] = getTrialsSU(data.stimEpocs,data.MUspikeT{cell_n}, cond, stim_duration);
        end

        spikes=Spike_Timing(:);
        spikes=spikes(spikes>latency & spikes<stim_duration+latency);
        R(cell_n,cond,rep) = length(spikes)/(stim_duration*numtrials);
        F1(cell_n,cond,rep) = spikeFFT(spikes,tf(rep));
        F2(cell_n,cond,rep) =spikeFFT(spikes,2*tf(rep));
        
        %%% rasters
        figure(spontfig);
        subplot(2,1,rep)
        hold on; set(gca, 'yDir','reverse');
        axis([0 plot_duration 0 numtrials+1]);
        plot ([Spike_Timing; Spike_Timing], [index-0.25;index+0.25], 'k', 'MarkerSize',4);
        set(gca,'XTickLabel',[]);     set(gca,'YTickLabel',[])
   
        if rep==1
            xlabel('spont')
            spont(cell_n) = R(cell_n,cond,rep);
        else
            xlabel('flicker')
        end
    end
    
    
    for rep =1:panels
        
        if SU
            rast_fig = figure('Name',sprintf('unit %d %d rep %d',channel_no,clust_no,rep));
        else
            rast_fig = figure('Name',sprintf('unit %d rep %d',cell_n));
        end
        
        for c =1:nrows*ncols;
            cond = (c-1)*panels+rep;
            if SU
                [Spike_Timing index numtrials] = getTrialsSU(stimEpocs{block},times, cond, stim_duration);
            else
                [Spike_Timing index numtrials] = getTrialsSU(data.stimEpocs,data.MUspikeT{cell_n}, cond, stim_duration);
            end
            %%% calculate total spikes
            spikes=Spike_Timing(:);
            spikes=spikes(spikes>latency & spikes<stim_duration+latency);
            R(cell_n,cond,rep) = length(spikes)/(stim_duration*numtrials);
            F1(cell_n,cond,rep) = spikeFFT(spikes,tf(rep))/(stim_duration*numtrials);
            F2(cell_n,cond,rep) =spikeFFT(spikes,2*tf(rep))/(stim_duration*numtrials) ;
            
            drawraster
      
   
            
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
            plot([0 max(hist_range)], [spont(cell_n) spont(cell_n)],'g')
            axis(axis_range);
            set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])
            
        end %orientation
        if SU
            saveas(rast_fig,fullfile(pname,sprintf('generic_rast_move%d%s_%d_%d',rep,Block_Name,channel_no,clust_no)),'fig');
            saveas(hist_fig,fullfile(pname,sprintf('generic_hist_move%d%s_%d_%d',rep,Block_Name,channel_no,clust_no)),'fig');
        end       
    end  %%% panel
    
        function drawraster
            figure(rast_fig);
            subplot(nrows,ncols,c)
            hold on; set(gca, 'yDir','reverse');
            axis([0 plot_duration 0 numtrials+1]);
            plot ([Spike_Timing; Spike_Timing], [index-0.25;index+0.25], 'k', 'MarkerSize',4);
            set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])
end
    
