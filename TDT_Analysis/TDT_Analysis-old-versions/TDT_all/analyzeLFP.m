% Matlab codes for reading from TTank

clear all;
pack
close all
Tank_Name='120205'
Block_Name='grating freq 1'
TTX = openTTX(Tank_Name,Block_Name); % to initialize

Event_Name_Snip='Snip'
Event_Name_Wave='PDec'
Event_Name_LFP = 'Lfpx'
Sample_Interval=0.04096 % 24414.0625Hz
Sample_Number_Snip=64
Dec_Factor=32; % 
plot_duration=8; %in second
hist_range=[0:0.1:9];
axis_range=[0 plot_duration 0 4];

tetrode_linear=1;

if tetrode_linear
    ch_map = [14 8 10 4 13 7 9 3 11 1 15 5 12 2 16 6];
else
    ch_map = 1:16;
end


figure
for ch=1:16
    %N = invoke(TTX, 'ReadEventsV', 5000, Event_Name_Snip, channel_no, 0,  Select_Duration(1), Select_Duration(2), 'ALL')
%    invoke(TTX,'CreateEpocIndexing');
% MyEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', 0,0, 1000);
% Select_Duration(1) = MyEpocs(2,1); % to exclude events before the first trigger      
Select_Duration(1)=0   
N = invoke(TTX, 'ReadEventsV', 10000, Event_Name_LFP, ch_map(ch), 0, Select_Duration(1),0,'ALL')
    W = invoke(TTX, 'ParseEvV', 0, N);
    event_times_all = invoke(TTX, 'ParseEvInfoV', 0, N, 6);  %   3  = event codes
    subplot(4,4,ch);
    plot( W(:,1:min(50,size(W,2))));
    axis([10 32 -10^(-4) 10^(-4)])
    axis off
    noise(ch) = std(W(1,:))/(10^-6);
    sprintf('channel %d : noise = %f',ch,noise(ch))
    phi(ch,:) = W(:);
end

CSD = phi(1:14,:) + phi(3:16,:) - 2*phi(2:15,:);
mean(std(CSD))
std(mean(CSD))

figure
for i = 1:12;
    subplot(6,2,i);
    %figure
    plot(CSD(i,:));
    axis([0 size(CSD,2) -10^-4 10^-4])
    %axis off
end
figure
plot(sum(CSD));
hold on
plot(std(CSD),'r');
std(sum(CSD))


plot_all = 1;
plot_chan =0;
if plot_all
    cluster_fig = figure;
    snip_fig = figure;
    hist_fig = figure;
end


max_events=20000;

idx_all = zeros(16,max_events);
event_times_all = zeros(16,max_events);
for ch=1:16
    

    ch_title = sprintf('channel %d',ch)
       invoke(TTX,'CreateEpocIndexing');
    MyEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', 0,0, 1000);
    Select_Duration(1) = MyEpocs(2,1); % to exclude events before the first trigger      
   % Select_Duration(1)=1300   
    N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_Snip, ch_map(ch), 0, Select_Duration(1),0,'ALL')
    W = invoke(TTX, 'ParseEvV', 0, N);
    event_times_all(ch,1:N) = invoke(TTX, 'ParseEvInfoV', 0, N, 6);  %   3  = event codes
    % figure
    % title(ch_title);
    % plot( W(:,1:50));
    % axis([10 32 -10^(-4) 10^(-4)])
    % axis off
    noise(ch) = std(W(1,:))/(10^-6);
    sprintf('channel %d : noise = %f',ch,noise(ch))

    X=W';

    X=double(X(:,13:42));
    [coeff score latent] = princomp(X(:,1:20));
    mn = mean(X);


    % figure
    % plot(latent);
    % title(ch_title);

    % figure
    % plot(score(:,1),score(:,2),'.');
    % figure
    % plot(score(:,1),score(:,3),'.');
    % figure
    % plot(score(:,2),score(:,3),'.');
    % figure
    % plot3(score(:,1),score(:,2),score(:,3),'.');

    %figure
    mn = mean(X);
    n_clust = 2;
% 100305_2
    if  (ch==14) | (ch==16)
        n_clust=2;
    end
    if (ch==9) | ch==8
        n_clust=3;
    end
    if ch ==13
        n_clust=4;
    end
% % 0930_4
%     if ch==7 | ch==8
%         n_clust=3;
%     end

% % 0930_6
% if ch==6 | ch==15
%     n_clust=3;
% end
    n_pca=3;
    if n_clust>1
        [idx c_mean] = kmeans(score(:,1:n_pca),n_clust);
    else
        idx = ones(size(X,1),1);
    end
    idx_all(ch,1:size(idx))=idx;

    %[idx] = clusterdata(score(:,1:n_pca),'linkage','ward','maxclust',3);
    %n_clust = max(idx)
    % for c = 1:n_clust
    %     plot(coeff(:,1:n_pca)*c_mean(c,:)'+mn');
    %     hold on
    % end

    group1 = X(find(idx==1),:);
    group2 = X(find(idx==2),:);
    group3 = X(find(idx==3),:);
    group4 = X(find(idx==4),:);

    

    clust1_size(ch) =size(group1,1);
    clust1_noise(ch) = std(min(group1'));
    clust2_size(ch) = size(group2,1);
    clust2_noise(ch) = std(min(group2'));

    skew1(ch) = skewness(min(group1'));
    skew2(ch) = skewness(min(group2'));

    hist_int = 10^-5 * [-15:0.2:0];

    if plot_all
        figure(hist_fig); 
        subplot(4,4,ch);
        hist_data = [hist(min(group1'),hist_int) ; hist(min(group2'),hist_int) ; ...
                hist(min(group3'),hist_int) ; hist(min(group4'),hist_int)];
        plot(hist_int*10^6,hist_data(1:n_clust,:)','linewidth',1.5);

        %hist(min(X'),hist_int); hold on
        set(gca,'XTickLabel',[])
       set(gca,'YTickLabel',[])
    else
        figure
        subplot(1,3,1);
        hist_data = [hist(min(group1'),hist_int) ; hist(min(group2'),hist_int) ; ...
                hist(min(group3'),hist_int) ; hist(min(group4'),hist_int)];
        plot(hist_int*10^6,hist_data(1:n_clust,:)','linewidth',1.5);
        axis([-150 0 0 max(max(hist_data))])
        xlabel('uVolt');
        if n_clust==2
            legend('1','2','Location','NorthWest');
        else
            legend('1','2','3','Location','NorthWest');
        end
        %hist(min(X'),hist_int); hold on
%         set(gca,'XTickLabel',[])
%         set(gca,'YTickLabel',[])
    end

    if plot_all
        figure(snip_fig); 
        subplot(4,4,ch);
    else
      
        subplot(1,3,2);
    end

    

    used = zeros(size(idx));
    %used(ceil(rand(100,1)*size(idx,1)))=1;
   used(1:(min(size(idx),75)))=1;
    if plot_all
        s_max = min(75,size(idx,1));
    else
        s_max = min(100,size(idx,1));
    end
    group1 = X(find((idx==1) & used),:)*10^6;
    group2 = X(find((idx==2) & used),:)*10^6;
    group3 = X(find((idx==3) & used),:)*10^6;
    group4 = X(find((idx==4) & used),:)*10^6;
    
    t = Sample_Interval*(0:size(group1,2)-1);
    plot(t, group1(1:min(100,size(group1,1)),:)','b','LineWidth',0.1);
    hold on
    plot(t,group2(1:min(100,size(group2,1)),:)','g','LineWidth',0.1);
    if n_clust>2
        plot(t,group3(1:min(100,size(group3,1)),:)','r','LineWidth',0.1);
    end
    if n_clust>3
        plot(t,group4(1:min(100,size(group4,1)),:)','k','LineWidth',0.1);
    end    
    axis([0 30*Sample_Interval -150 100])
      if plot_all
          set(gca,'XTickLabel',[])
          set(gca,'YTickLabel',[])
      end
      if plot_chan
          ylabel('uVolt');
          xlabel('msec');
      end


    used = zeros(size(idx));
    %used(ceil(rand(1000,1)*size(idx,1)))=1;
    used(1:(min(size(idx),500)))=1;
       s_max = min(500,size(idx,1));
    if plot_all
        figure(cluster_fig);
        subplot(4,4,ch);
    else
        subplot(1,3,3);
    end
    plot(score(find(idx==1 & used),1),score(find(idx==1 & used),2),'.','MarkerSize',1);
    hold on
    plot(score(find(idx==2 & used),1),score(find(idx==2 & used),2),'g.','MarkerSize',1);
    if n_clust>2
        plot(score(find(idx==3 & used),1),score(find(idx==3 & used),2),'r.','MarkerSize',1);
    end
    if n_clust>3
        plot(score(find(idx==4 & used),1),score(find(idx==4 & used),2),'k.','MarkerSize',1);
    end
    axis((10^-4)*[-2 2 -2 2]);
    set(gca,'XTickLabel',[])
    set(gca,'YTickLabel',[])    
    
    if plot_chan
        set(gcf,'Position',[180   435   757   235]);
            title(ch_title);

    end
    
    %%% calculate ISI
    
    ISI_hist_range = (0:2:100)*(10^-3);
    ISI_hist = zeros(n_clust, size(ISI_hist_range,2));
    for cl = 1:n_clust
        spike_times=squeeze(event_times_all(ch,find(idx==cl)));
        ISI=spike_times(2:size(spike_times,2)) - spike_times(1:size(spike_times,2)-1);
        ISI_hist(cl,:) = hist(ISI,ISI_hist_range);
        min(ISI)
    end
    
    figure
    n_pts = 50;
    plot(ISI_hist_range(1:n_pts),ISI_hist(:,1:n_pts));
   % axis ([0 ISI_hist_range(n_pts) 0 max(max(ISI_hist(:,1:n_pts)))])
    hold on
    spike_times = squeeze(event_times_all(ch,find(idx~=0)));
    ISI = spike_times(2:size(spike_times,2)) - spike_times(1:size(spike_times,2)-1);
   min(ISI)
    ISI_hist(cl,:) = hist(ISI,ISI_hist_range);
       plot(ISI_hist_range(1:n_pts),ISI_hist(:,1:n_pts));
       
    
    
%         figure
%         title(ch_title);
%         plot3(score(find(idx==1),1),score(find(idx==1),2),score(find(idx==1),3),'.');
%         hold on
%         plot3(score(find(idx==2),1),score(find(idx==2),2),score(find(idx==2),3),'g.');
%         hold on
%         plot3(score(find(idx==3),1),score(find(idx==3),2),score(find(idx==3),3),'r.');
%         title(ch_title);
  
end

% 
% 
% fname='pca.fet.1'
% dlmwrite(fname,n_pca);
% dlmwrite(fname,score(:,1:n_pca),'-append','delimiter',' ');
% 
% cl_name = 'pca.clu.1';
% idx = dlmread(cl_name);
% n_clust = idx(1);
% idx = idx(2:size(idx,1));
% 
% figure
% group1 = X(find(idx==1),:);
% group2 = X(find(idx==2),:);
% group3 = X(find(idx==3),:);
% 
% subplot(4,1,1);
% plot(group1(1:100,:)','b');
% subplot(4,1,2);
% plot(group2(1:100,:)','r');
% subplot(4,1,3);
% plot(group3(1:100,:)','g');
% 
% 
% figure
% plot3(score(find(idx==1),1),score(find(idx==1),2),score(find(idx==1),3),'.');
% hold on
% plot3(score(find(idx==2),1),score(find(idx==2),2),score(find(idx==2),3),'g.');
% hold on
% plot3(score(find(idx==3),1),score(find(idx==3),2),score(find(idx==3),3),'r.');
% 
% figure
% hold on
% plot(X(find(idx==1),:)','b');
% plot(X(find(idx==2),:)','r');
% plot(X(find(idx==3),:)','g');
% plot(X(find(idx==4),:)','k');
% 
% 
% TimeStamp = invoke(TTX, 'ParseEvInfoV', 0, N, 6);
% % save (Save_File_Snip, 'TimeStamp', 'W');
% 
% ISI =  diff(TimeStamp);
% figure; Hist(ISI, [0:0.0025:0.5]); axis ([-0.02 0.2 0 100]);
% 
% % N = invoke(TTX, 'ReadEventsV', 10000, 'xTrg', 0, 0, select_duration(1), select_duration(2), 'ALL')
% % xTrglist = invoke(TTX, 'ParseEvV', 0, N);
% % xTrg_TS = invoke(TTX, 'ParseEvInfoV', 0, N, 6);
% 
% invoke(TTX, 'CloseTank')
% invoke(TTX, 'ReleaseServer');
% 
% % To get other information about the record events returned call
% % ParseEvInfoV.  This call has the same two parameters as ParseEvV
% % with one more param to indicate which bit of information you
% % want returned.  The following are valid values for the 3rd 
% % parameter:
% %   1  = Amount of waveform data in bytes
% %   2  = Record Type (see TCommon.h)
% %   3  = Event Code Value
% %   4  = Channel No.
% %   5  = Sorting No.
% %   6  = Time Stamp
% %   7  = Scalar Value (only valid if no data is attached)
% %   8  = Data format code (see TCommon.h)
% %   9  = Data sample rate in Hz. (not value unless data is attached)
% %   10 = Not used returns 0.
% %   0  = Returns all values above