%function cluster_linear_multexp
% clusters data from linear multi-electrode recording
% reads in snippets, aligns them, calculates principal components
% sends PC to Klustakwik (dos program, from Buzsaki)
% then displays various diagnostics (clusters, amplitudes, waveforms)

% cmn 06-06

clear all;
pack
close all
Tank_Name='061107_wt_linear'
nblocks =1;

Block_Name={'drift2'  'wn2' 'sb2' 'cp_tf2' 'flash_grat1' 'bars16d2' 'bars16d2b'}
nblocks = size(Block_Name,2)

TTX = openTTX(Tank_Name,char(Block_Name(1))); % to initialize

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
thresh = -10^-6 * ones(16,1) * 30;
thresh(10:15) = -40* 10^-6;

%%% create figuresw
cluster_fig = figure;
snip_fig = figure;
hist_fig = figure;
cluster_fig2 = figure;



%%% set time limits
max_events=30000;
start_time = 60*0;   %% in seconds
max_time = 60*30;   

idx_all = zeros(16,max_events);
event_times_all = zeros(16,max_events);
for ch=1:16
   clear X
   pack
   
    ch_title = sprintf('channel %d',ch)
    N=0;
    Nblock=0;
    for block= 1:nblocks;
 
       TTX = openTTX(Tank_Name,char(Block_Name(block))); % to initialize
        invoke(TTX,'CreateEpocIndexing');
        MyEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', 0,0, 1000);
        Select_Duration(1) = max(MyEpocs(1,1),start_time); % to exclude events before the first trigger      

        Nblock(block) = invoke(TTX, 'ReadEventsV', max_events, Event_Name_Snip, ch_map(ch), 0, Select_Duration(1),max_time,'ALL')
        W = invoke(TTX, 'ParseEvV', 0, Nblock(block));
        event_times_all(ch,N+(1:Nblock(block))) = invoke(TTX, 'ParseEvInfoV', 0, Nblock(block), 6) +10^5 *(block-1);      %   3  = event codes

        noise(ch) = std(W(1,:))/(10^-6);
        sprintf('channel %d : noise = %f',ch,noise(ch))
        X(N+(1:Nblock(block)),1:size(W,1))=W';
        blockID(ch,N+(1:Nblock(block)))=block;
        N=N+Nblock(block);
    end
    %%% align spikes by bottom of the trough, and eliminate those that
    %%% don't meet threshold
    clear W
    
    [Y I] = min(X(:,8:14),[],2);
    used = find((Y<thresh(ch)) );
    cutoff = max(Y)
    X = X(used,:);
    N = size(used)
    I = I(used);
    tic
    Xshift = zeros(N,22);
    for i = 1:N
        Xshift(i,:) = X(i,I(i) : (I(i)+21));
    end
    toc
    X = Xshift;
    clear Xshift

   %%% calculate independent components and send to KlustaKwik
   [coeff score latent] = princomp(X);
        n_pca=4;
    [s coeff u] = fastica(X'*10^9,'numOfIC',n_pca,'lastEig',n_pca,'g','tanh','stabilization','on');
       score = s';

    
    
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
        mean_wvform(:,ch, c) = mean(group(1:size(members,1),:,c));
    end
   
%%% plot histograme of spike amplitude (i.e. minimum value)
    hist_int = 10^-5 * [-15:0.2:0];
    figure(hist_fig); 
    subplot(4,4,ch);        
      for cl = 1:n_clust
          hist_data(cl,:) = hist(min(squeeze(group(1:clust_size(cl),6:10,cl)),[],2),hist_int);
      end
    plot(hist_int*10^6,hist_data(1:(min(n_clust,7)),:)','linewidth',1.5);
    axis([min(hist_int)*10^6 max(hist_int)*10^6 0 200 ]);
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

     np =5000;    %% select np points to display
    N = size(score,1);
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

         np =5000;    %% select np points to display
        N = size(score,1);
        used=zeros(size(idx));
        used(ceil(rand(np,1)*N))=1;

        for cl = 1:n_clust
            plot(10^6*score(find(idx==cl & used),3),10^6*score(find(idx==cl & used),4),colorstyle((cl*2-1):cl*2),'MarkerSize',1);
            hold on
        end
        
    set(gca,'XTickLabel',[])
    set(gca,'YTickLabel',[]) 
    end
    
    
  for c= 1:n_clust
      c_score = score(find(idx==c),:);
      mn = mean(c_score,1);
      score_cov = cov(c_score);
      nonc_score = score(find(idx~=c),:);
      cov_inv = inv(cov(c_score));
      
      Dnon = zeros(1,size(nonc_score,1));
      Dc = zeros(1,size(c_score,1));
      var_inv = 1./var(c_score);

%       for i = 1:size(c_score,1)
%           Dc(i) = (c_score(i,:)-mn)*cov_inv*(c_score(i,:)-mn)';
%       end
%       for i = 1:size(nonc_score,1)
%           Dnon(i) = (nonc_score(i,:)-mn)*cov_inv*(nonc_score(i,:)-mn)';
%       end


       for i = 1:size(c_score,1)
           Dc(i) = (var_inv.*(c_score(i,:)-mn))*(c_score(i,:)-mn)';
       end
       for i = 1:size(nonc_score,1)
           Dnon(i) = (var_inv.*(nonc_score(i,:)-mn))*(nonc_score(i,:)-mn)';
       end
      
        L(ch,c) = sum(1-chi2cdf(Dnon,size(score,2)));
      Lratio(ch,c) = sum(1-chi2cdf(Dnon,size(score,2)))/size(c_score,1);
      csize(ch,c) = size(c_score,1);
  end
L
   Lratio
  csize 
    

  for c = 1:n_clust %%% waveform parameters
     c
%      figure
%      plot(1:size(mean_wvform,1),mean_wvform (:,(tet-1)*4+1 : (tet-1)*4+4,c)');
    trig_chan = ch; %%% legacy from tetrodes
  
    
    
       [trough(ch,c) trough_t] = min(mean_wvform (:,ch,c));
      [peak(ch,c) peak_t] = max(mean_wvform ((trough_t+1):size(mean_wvform,1),trig_chan,c));
    
          if trough_t+peak_t<size(mean_wvform,1)
           y1 = mean_wvform(trough_t+peak_t-1,trig_chan,c);
          y2 = mean_wvform(trough_t+peak_t,trig_chan,c);
          y3 = mean_wvform(trough_t+peak_t+1,trig_chan,c);
          %%% fit parabola through three points to find trough
          peak_t = peak_t - 0.5*((y3-y1)/(y3+y1-2*y2));
      end   
      
      if trough_t>1
          y1 = mean_wvform(trough_t-1,trig_chan,c);
          y2 = mean_wvform(trough_t,trig_chan,c);
          y3 = mean_wvform(trough_t+1,trig_chan,c);
          %%% fit parabola through three points to find trough
          trough_t =  - 0.5*((y3-y1)/(y3+y1-2*y2));
      end
      
     
          
      halfheight = trough(ch,c)/2;
    trough_points = find(mean_wvform(:,trig_chan,c)<halfheight);
    half_minus = min(trough_points);
    if half_minus>1
        half_minus = half_minus - (halfheight - mean_wvform(half_minus,trig_chan,c))/(mean_wvform(half_minus-1,trig_chan,c) - mean_wvform(half_minus,trig_chan,c));
    end
    half_plus = max(trough_points);
    if half_plus<size(mean_wvform,1)
        halfplus = half_plus + (halfheight - mean_wvform(half_plus,trig_chan,c))/(mean_wvform(half_plus+1,trig_chan,c) - mean_wvform(half_plus,trig_chan,c));
    end
    half_width(ch,c) = half_plus - half_minus;
  peak_trough(ch,c) = peak_t-trough_t;
  
  end

  clear s score score_cov c_score nonc_score
end


%%% save results if so desired
%%% saves all variables, plus figures
Block_Name
resp=input('save results y/n','s')
if resp=='y'
      bname = input('filemaes : ','s')
    output_path=uigetdir('','data folder')
    fname = fullfile(output_path,sprintf('clustered%s_%s',Tank_Name,bname));
    save(fname);
    fname = fullfile(output_path,sprintf('hist%s_%s',Tank_Name,bname));
    saveas(hist_fig,fname,'fig');
     fname = fullfile(output_path,sprintf('snip%s_%s',Tank_Name,bname));
    saveas(snip_fig,fname,'fig');   
    fname = fullfile(output_path,sprintf('clust%s_%s',Tank_Name,bname));
    saveas(cluster_fig,fname,'fig');
        fname = fullfile(output_path,sprintf('clust2%s_%s',Tank_Name,bname));
    saveas(cluster_fig2,fname,'fig');
end

