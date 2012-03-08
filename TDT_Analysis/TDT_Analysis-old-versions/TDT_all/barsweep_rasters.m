%function barsweep16d_cluster
% Matlab codes for reading from TTank for sweeping bars in 8 orientations
% plots histgrams and rasters and fits data to a gaussian peak
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03



%%% read in cluster data, then connect to the tank and read the block
clear all
[fname, pname] = uigetfile('*.mat','cluster data');
if fname~=0
    load(fullfile(pname,fname));
    if exist('Block_Name')
          block =6;
         Block_Name = char(Block_Name(block)) % to initialize
    else
        block=2;
        Block_Name = char(Block_Name(block))
    end
    use_clusters=1;
else
    use_clusters=0;
    Tank_Name='012507_wt_linear'
    Block_Name={'bars16d1' }
end

tetrode_linear=0;

if tetrode_linear
    ch_map = [14 8 10 4 13 7 9 3 11 1 15 5 12 2 16 6];
else
    ch_map = 1:16;
end


[afname, apname] = uigetfile('*.mat','analysis data');
if afname~=0
  afile = fullfile(apname,afname);
  load(afile); 
  use_afile=1;
else
    use_afile=0;
    %%% select which units to analyze (channel, cluster number)
 
    cells = [1 1; 1 2; 5 1; 5 2; 5 4; 5 6; 5 7; 9 4; 9 5;  9 3; 9 7; 13 4; 13 3; 13 7]; %%%011707_allstim1
    cells = [8 1; 8 2; 10 1; 11 3; 12 3; 12 2; 12 6; 13 4; 13 2; 13 3; 14 3; 14 5; 15 2; 15 3; 15 4]
    cells = [7 3]
end

if ~use_clusters
    cells = [1 1; 2 1; 3 1; 4 1; 5 1 ; 6 1; 7 1; 8 1; 9 1; 10 1; 11 1; 12 1; 13 1; 14 1; 15 1; 16 1];
end    


event_times_all = event_times_all-(block-1)*10^5;
times = event_times_all(event_times_all>0 & event_times_all<10^5);

   linecolor = [0 0 1; 0 1 0 ; 1 0 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0; .25 0 0.5; 0 0.5 0 ; 0.5 .25 0; 0.5 0 1; 0 0.5 0.5];

figure
for cell_n = 1:size(cells,1)

    cell_n
    channel_no = cells(cell_n,1)
    clust_no = cells(cell_n,2)

  channel_times =squeeze(event_times_all(channel_no,:));

spike_n =find(channel_times>0 & channel_times<60 & idx_all(channel_no,:)==clust_no);
    t = channel_times(spike_n);
    
    c = ones(size(t));
    c(:)= cell_n;
    hold on
    for i = 1:length(t)
        plot([t(i) t(i)],[(cell_n-0.35) (cell_n+0.35)],'Color',linecolor(cell_n,:),'LineWidth',1);
    end 
end    
