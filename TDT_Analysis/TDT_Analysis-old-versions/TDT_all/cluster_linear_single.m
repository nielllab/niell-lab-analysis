%function cluster_linear
% clusters data from linear multi-electrode recording
% reads in snippets, aligns them, calculates principal components
% sends PC to Klustakwik (dos program, from Buzsaki)
% then displays various diagnostics (clusters, amplitudes, waveforms)

% cmn 06-06

clear all;
pack
close all
Tank_Name='091806_wt'
Block_Name='spot10deg'
TTX = openTTX(Tank_Name,Block_Name); % to initialize

Event_Name_Snip='Snip'
Sample_Interval=0.04096 % 24414.0625Hz

%%% remap site numbers if using a connector meant for tetrode configuration
%%% not used since 9/05, but could be used for other remappings
tetrode_linear=0;
if tetrode_linear
    ch_map = [14 8 10 4 13 7 9 3 11 1 15 5 12 2 16 6];
else
    ch_map = 1:16;
end

%%% set a software threshold (minimum amplitude) for snippets
%%% can be different for each channel
thresh = -10^-6 * ones(16,1) * 10;
% thresh(10) = -35* 10^-6;
% thresh(11) = -35* 10^-6;
% thresh(12) = -35* 10^-6;
% thresh(16) = -35* 10^-6;

%%% create figuresw
cluster_fig = figure;
snip_fig = figure;
hist_fig = figure;
cluster_fig2 = figure;



%%% set time limits
max_events=50000;
start_time = 60*0;   %% in seconds
max_time = 60*30;   

idx_all = zeros(16,max_events);
event_times_all = zeros(16,max_events);
for ch=1:1
    ch_title = sprintf('channel %d',ch)
       invoke(TTX,'CreateEpocIndexing');
    MyEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', 0,0, 1000);
    Select_Duration(1) = max(MyEpocs(1,1),start_time); % to exclude events before the first trigger      
 
    N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_Snip, ch_map(ch), 0, Select_Duration(1),max_time,'ALL')
    W = invoke(TTX, 'ParseEvV', 0, N);
    X=W';
    event_times_all(ch,1:N) = invoke(TTX, 'ParseEvInfoV', 0, N, 6);  %   3  = event codes

    noise(ch) = std(W(1,:))/(10^-6);
    sprintf('channel %d : noise = %f',ch,noise(ch))

    %%% align spikes by bottom of the trough, and eliminate those that
    %%% don't meet threshold
    
    [Y I] = max(X(:,8:12),[],2);
    used = find((Y>thresh(ch)) );
    cutoff = min(Y)
    X = X(used,:);
    N = size(used)
    I = I(used);
    I(:)=1;
    tic
    Xshift = zeros(N,26);
    for i = 1:N
        Xshift(i,:) = X(i,I(i) : (I(i)+25));
    end
    toc
    X = Xshift;

   %%% calculate principal components and send to KlustaKwik
   [coeff score latent] = princomp(X);
    n_pca=2;
    fname='pca.fet.1'
    dlmwrite(fname,n_pca);
    dlmwrite(fname,score(:,1:n_pca),'-append','delimiter',' ');
    dos('C:\KlustaKwik pca 1 -MaxClusters 8 -MaxPossibleClusters 8 -Verbose 0');
    cl_name = 'pca.clu.1';
    idx = dlmread(cl_name);

   %%% first cluster is outliers (noise), so get rid of these
   n_clust = idx(1)
    if n_clust>1
        n_clust=n_clust-1;
        idx = idx-1;
    end
    idx = idx(2:size(idx,1));
    idx_all(ch,used)=idx;
 
    %%% create an array with each cluster in a row, for easy graphing
    for c = 1:n_clust
        members = find(idx==c);
        clust_size(c)=size(members,1);
        group(1:size(members,1),1:size(X,2),c)=X(members,:);
    end
   
%%% plot histograme of spike amplitude (i.e. minimum value)
    hist_int = 10^-5 * [-10:0.2:10];
    figure(hist_fig); 
    subplot(4,4,ch);        
      for cl = 1:n_clust
          hist_data(cl,:) = hist(min(squeeze(group(1:clust_size(cl),10:25,cl)),[],2),hist_int);
      end
    plot(hist_int*10^6,hist_data(1:(min(n_clust,7)),:)','linewidth',1.5);
    axis([min(hist_int)*10^6 max(hist_int)*10^6 0 100 ]);
    set(gca,'XTickLabel',[])
    set(gca,'YTickLabel',[])



%%% plot waveforms
   figure(snip_fig); 
    subplot(4,4,ch);
    linecolor = 'bgrcmyk';
    t = (0:size(group,2)-1);
    np = 5; 
    for cl = 1:n_clust;
         plot(t, group(ceil(rand(np,1)*clust_size(cl)),:,cl)'*10^6,linecolor(cl),'LineWidth',0.1);
        hold on
    end

    axis([0 25 -100 75]) 
    set(gca,'XTickLabel',[])
    set(gca,'YTickLabel',[])



%%% cluster plot of 1st vs 2nd principal component

    figure(cluster_fig);
    subplot(4,4,ch);
    colorstyle = ['b.g.r.c.m.y.k.'];

     np =1000;    %% select np points to display
    N = min(np,size(score,1));
    used=zeros(size(idx));
    used(ceil(rand(np,1)*N))=1;
    
    for cl = 1:n_clust
        plot(10^6*score(find(idx==cl & used),1),10^6*score(find(idx==cl & used),2),colorstyle((cl*2-1):cl*2),'MarkerSize',1);
        hold on
    end
    set(gca,'XTickLabel',[])
    set(gca,'YTickLabel',[]) 
    
    %%% same for 3 vs 4
    if n_pca>2

        figure(cluster_fig2);
        subplot(4,4,ch);
        colorstyle = ['b.g.r.c.m.y.k.'];

         np =1000;    %% select np points to display
        N = min(np,size(score,1));
        used=zeros(size(idx));
        used(ceil(rand(np,1)*N))=1;

        for cl = 1:n_clust
            plot(10^6*score(find(idx==cl & used),3),10^6*score(find(idx==cl & used),4),colorstyle((cl*2-1):cl*2),'MarkerSize',1);
            hold on
        end
        
    set(gca,'XTickLabel',[])
    set(gca,'YTickLabel',[]) 
    end
    
    
%     %%% calculate ISI
%     
%     ISI_hist_range = (0:2:100)*(10^-3);
%     ISI_hist = zeros(n_clust, size(ISI_hist_range,2));
%     for cl = 1:n_clust
%         spike_times=squeeze(event_times_all(ch,find(idx==cl)));
%         ISI=spike_times(2:size(spike_times,2)) - spike_times(1:size(spike_times,2)-1);
%         ISI_hist(cl,:) = hist(ISI,ISI_hist_range);
%         min(ISI)
%     end
 
    
    
%     figure
%     n_pts = 50;
%     plot(ISI_hist_range(1:n_pts),ISI_hist(:,1:n_pts));
%    % axis ([0 ISI_hist_range(n_pts) 0 max(max(ISI_hist(:,1:n_pts)))])
%     hold on
%     spike_times = squeeze(event_times_all(ch,find(idx~=0)));
%     ISI = spike_times(2:size(spike_times,2)) - spike_times(1:size(spike_times,2)-1);
%    min(ISI)
%     ISI_hist(cl,:) = hist(ISI,ISI_hist_range);
%        plot(ISI_hist_range(1:n_pts),ISI_hist(:,1:n_pts));
       

  
end


%%% save results if so desired
%%% saves all variables, plus figures
Block_Name
resp=input('save results y/n','s')
if resp=='y'
    output_path=uigetdir('','data folder')
    fname = fullfile(output_path,sprintf('clustered%s_%s',Tank_Name,Block_Name));
    save(fname);
    fname = fullfile(output_path,sprintf('hist%s_%s',Tank_Name,Block_Name));
    saveas(hist_fig,fname,'fig');
     fname = fullfile(output_path,sprintf('snip%s_%s',Tank_Name,Block_Name));
    saveas(snip_fig,fname,'fig');   
    fname = fullfile(output_path,sprintf('clust%s_%s',Tank_Name,Block_Name));
    saveas(cluster_fig,fname,'fig');
end

