%function drift_analysis
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

stim_duration = input('duration [1] : ');
if stim_duration==[]
    stim_duration=1;
end
nrows = input('# orients [8] : ');

ncols = input('# sfs [6]  : ');

tf = input('temp freqs [2 8] : ');

latency = input('latency (0.05) : ');

latency
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
%for cell_n=21:26
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
        
        %%% rasters
        figure(spontfig);
        subplot(2,1,rep)
        hold on; set(gca, 'yDir','reverse');
        axis([0 plot_duration 0 numtrials+1]);
        plot ([Spike_Timing; Spike_Timing], [index-0.25;index+0.25], 'k', 'MarkerSize',4);
        set(gca,'XTickLabel',[]);     set(gca,'YTickLabel',[])
        
        if rep==1
            xlabel('spont')
            for f = 0:2
                drift(cell_n).spont(f+1) = spikeFFT(spikes,tf(1)*f)/(stim_duration*numtrials);
            end
        else
            for f = 0:2  %%% fix this when adding more tf for flicker
                drift(cell_n).flicker(f+1) = spikeFFT(spikes,tf(1)*f)/(stim_duration*numtrials);
            end
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
            sf_ind = mod(c-1,ncols)+1;
            orient_ind= ceil(c/ncols);
            if SU
                [Spike_Timing index numtrials] = getTrialsSU(stimEpocs{block},times, cond, stim_duration);
            else
                [Spike_Timing index numtrials] = getTrialsSU(data.stimEpocs,data.MUspikeT{cell_n}, cond, stim_duration);
            end
            %%% calculate total spikes
            spikes=Spike_Timing(:);
            spikes=spikes(spikes>latency & spikes<stim_duration+latency);
            for f = 0:2
                drift(cell_n).R(orient_ind,sf_ind,rep,f+1) = spikeFFT(spikes,tf(rep)*f)/(stim_duration*numtrials);
            end
            
            figure(rast_fig);
            subplot(nrows,ncols,c)
            hold on; set(gca, 'yDir','reverse');
            axis([0 plot_duration 0 numtrials+1]);
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
            plot([0 max(hist_range)], [drift(cell_n).spont(1) drift(cell_n).spont(1)],'g')
            axis(axis_range);
            set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])
            
        end %orientation
        if SU
            saveas(rast_fig,fullfile(pname,sprintf('drift_rast_move%d%s_%d_%d',rep,Block_Name,channel_no,clust_no)),'fig');
            if rep==2
                saveas(hist_fig,fullfile(pname,sprintf('drift_hist_move%d%s_%d_%d',rep,Block_Name,channel_no,clust_no)),'fig');
            end
        end
        if rep==1
            colorplot=figure
        else
            figure(colorplot)
        end
        for f = 1:3
            subplot(2,3,f+(rep-1)*3)
            imagesc(squeeze(abs(drift(cell_n).R(:,:,rep,f))));
            colorbar
        end
        
    end  %%% panel
    
    %wfig= figure
    thetafig=figure('Name',sprintf('unit %d %d rep %d',channel_no,clust_no,rep));
    color={'b','r'};
    drift(cell_n).orient_tune = zeros(3,length(tf),nrows);
    drift(cell_n).sf_tune=zeros(3,length(tf),ncols+1);
    
    for f=1:3
        for rep = 1:length(tf)
           
            [u s v] = svd(abs(squeeze(drift(cell_n).R(:,:,rep,f)))-abs(drift(cell_n).spont(f)))
           %[u s v] = svd(abs(squeeze(drift(cell_n).R(:,:,rep,f))))
          
           orient_tune = u(:,1);
           
           sf_tune = v(:,1);
                      if sum(orient_tune)<0 & sum(sf_tune)<0;
                orient_tune=-1 * orient_tune;
                sf_tune=-1 * sf_tune;
                      end 
            [max_o]= max(orient_tune);
            max_sf = max(sf_tune);
                sf_tune = sf_tune*s(1,1)*max_o;
                orient_tune=orient_tune*s(1,1)*max_sf;
                
                
           sf_tune(2:end+1)=sf_tune;
            sf_tune(1) = abs(drift(cell_n).flicker(f));  %%% fix this when adding more tfs for flicker
            
            drift(cell_n).orient_tune(rep,f,:)=orient_tune';
            drift(cell_n).sf_tune(rep,f,:) = sf_tune;
 
            figure(thetafig)
            subplot(2,3,f)
            hold on
            plot(orient_tune,color{rep})
            xlabel('theta')
            lim = ylim;
            if rep==2
                ylim([min(0,lim(1)) lim(2)])
            end 
            xlim([1 length(orient_tune)])
            
            %figure(wfig)
            subplot(2,3,f+3)
            hold on
            plot(sf_tune,color{rep})
            xlabel('SF')
            xlim([1 length(sf_tune)])
            if rep ==2
                lim = ylim;
            ylim([min(0,lim(1)) lim(2)])
            end
        end
        
    end
end

save(afile,'drift','-append');
