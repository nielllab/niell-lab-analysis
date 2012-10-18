%function timecourse_analysis
% Matlab codes for reading from TTank for sweeping bars in 8 orientations
% plots histgrams and rasters and fits data to a gaussian peak
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03



%%% read in cluster data, then connect to the tank and read the block
clear all
SU = input('multiunit (0) or single-unit (1) : ');

cells =1;
if SU
    [fname, pname] = uigetfile('*.mat','cluster data');
    load(fullfile(pname,fname));
    for i =1:length(Block_Name);
        sprintf('%d : %s ',i,Block_Name{i})
    end
    block = input('which block to analyze ? ');
    Block_Name = Block_Name{block}
    [afname, apname] = uigetfile('*.mat','analysis data');
    
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

binsize = input('bin size (sec) : ');

if SU
    cell_range = 1:size(cells,1)
else
    cell_range=1:4:nchan;
end
for cell_n = cell_range;
    % for cell_n=2:2
    % for cell_n=9:9
    cell_n
    if SU
        channel_no = cells(cell_n,1)
        clust_no = cells(cell_n,2)
        channel_times =spikeT{cell_n} - (block-1)*10^5;
        use = find(channel_times>0 & channel_times<10^5);
        wvall = wave_all{ceil(channel_no/4)};
        wvclust = wvall(find(idx_all{channel_no}==clust_no),:,:);
        
        times = channel_times(use);
        wv = wvclust(use,:,:);
        amps =squeeze(min(wv(:,5:10,:),[],2));
        figure
        hist(amps,[-200:1:0]*10^-6)
        if SU
            title(sprintf('ch %d cl %d',channel_no,clust_no));
        else
            title(sprintf('channel %d',cell_n));
        end
        
        figure
        plot(times/60,amps,'.');
        if SU
            title(sprintf('ch %d cl %d',channel_no,clust_no));
        else
            title(sprintf('channel %d',cell_n));
        end
        
    else
        times = data.MUspikeT{cell_n};
    end
    
    
    
    
    histbins = (0:binsize:max(times)-binsize/2) + binsize/2;
    h= hist(times,histbins);
    figure
    plot(histbins/60,h);
    if SU
        title(sprintf('ch %d cl %d',channel_no,clust_no));
    else
        title(sprintf('channel %d',cell_n));
    end
    axis([0 max(histbins)/60 0 max(h)*1.1])
    
end

