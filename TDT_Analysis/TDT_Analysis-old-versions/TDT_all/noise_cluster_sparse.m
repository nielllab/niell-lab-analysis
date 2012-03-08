%function noise_cluster

% Matlab codes for reading from TTank for movie data 
% Calculates RFs from spike triggered average, and spike triggered covariance
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03


%%% read in cluster data, then connect to the tank and read the block
clear all
pack
[fname, pname] = uigetfile('*.mat','cluster data');
oldpname = pname;
load(fullfile(pname,fname));
pname = oldpname;
Block_Name
if nblocks>1
      block =4;
     Block_Name = char(Block_Name(block)) % to initialize
else
    block=2;
end
%load wn012cpd10Hz075sig5min_0417a   %%% movie frames
%load noisewhite012cpd3hz2std3min
%load noisewhite012cpd5Hz3min_0308

%load wn_movie_5hz_30hz_012cpd_075sig_5min_092106
%load wn012cpd75hz30hz075sig5min_101206
%load wn016cpd30hz75hz075sig3min_011607
% load wn08cpd5hz075sig30hz5min_012307
% load wn012alpha1_5hz_contrast10sec5min021407   %%% the usual
% n_frames = 8999;
% load 012flat_nocmod_2to8Hz_05sig_10min
% n_frames = 17999;
%load sparse5deg5percent2x052608
 %load sparse5deg05percent5mindouble
 %load sparse5deg10percentdouble
%load sparse5deg20percent
load sparse10deg05perdouble
n_frames=17999
sta_size = 4
filt = fspecial('gaussian',10,sta_size);
filt = filt/sum(sum(filt));
%load wn012cpd5hz30hz_alph1_off3_5min_021207
%load alpha_stepcontrast
%load sf_movie
clear amp width x0 baseline
TTX = openTTX(Tank_Name,Block_Name); % to initialize

Event_Name_Snip='Snip'
Event_Name_Wave='PDec'
Sample_Interval=0.04096; % 24414.0625Hz
Sample_Number_Snip=64;
Dec_Factor=32; % 
plot_duration=36; %in second
hist_int = 0.1;
hist_range=[0:hist_int:plot_duration];
axis_range=[0 plot_duration 0 20];
max_events=100000;
 xrange =25:45;
yrange = 30:50;   
dx = size(xrange,2);
dy=size(yrange,2);

moviedata = double(moviedata);
mov_var = var(squeeze(moviedata(50,50,:)))
movavg = mean(moviedata,3);

figure
imagesc(movavg-128, [-32/1.414 32/1.414])



[fname, pname] = uigetfile('*.mat','analysis data');

  afile = fullfile(pname,fname);
  load(afile); 
  use_afile=1;

%%% set time based on first and last timepoints in clustered data
invoke(TTX,'CreateEpocIndexing')

event_times_all = event_times_all-(block-1)*10^5;
times = event_times_all(event_times_all>0 & event_times_all<10^5);
 time1 = min(times);
 time2= max(times);


MyEpocs = invoke(TTX, 'GetEpocsV', 'fTrg',time1, time2, round(max_time*30));

printfig = input('print figures ? ');

for cell_n = 1:size(cells,1)

    channel_no = cells(cell_n,1)
    clust_no = cells(cell_n,2)

    n_spikes= zeros(n_frames,1);
    
   
  channel_times =squeeze(event_times_all(channel_no,:));
  times = channel_times(channel_times>0 & channel_times<10^5);
  time1 = min(times)
  time2=max(times)
    
    sta0 = zeros(size(moviedata(:,:,1)));
    sta1 = sta0; sta2=sta0;

        clear Epocs trial_no Epocs_TS N Spike_TS index TS_xTrg Spike_Timing;
        invoke(TTX,'ResetFilters');
       N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_Snip, ch_map(channel_no), 0, ...
                    time1, time2, 'FILTERED')
        if (N==max_events)
            warning('max number of events aquired');
        end
        Spike_TS = invoke(TTX, 'ParseEvInfoV', 0, N, 6);  %   6  = Time Stamp
        n_reps = max(Spike_TS)/(.03*n_frames)
        
       
        %%% match up spike times for this condition with spike times from clustering
       %%% to identify cluster index for each spike     
%         channel_times =squeeze(event_times_all(channel_no,:));
%         spike_no=zeros(N,1);        
%         n=1;
%         for i = 1:N        
%             while channel_times(n)~=Spike_TS(i)
%                 n = n+1;
%             end
%             spike_no(i) = n;
%         end        
%         clust_idx = idx_all(channel_no,spike_no);

      channel_times =squeeze(event_times_all(channel_no,:));
        ntimes = size(event_times_all,2);
        spike_no=zeros(N,1);        
        n=1;
        for i = 1:N        
            while channel_times(n)~=Spike_TS(i) & n<ntimes
                n = n+1;
            end
            if channel_times(n)==Spike_TS(i)
                spike_no(i) = n;
            else    %% if not caught, start over at beginning for next search
                n=1;
                spike_no(i)=ntimes+1;
            end
                
        end        
       idx_all(channel_no,ntimes+1)=0;  %%%catchall for spikes not identified
       clust_idx = idx_all(channel_no,spike_no);
     
        %%% found out which frame each spike occured during
        
        Spike_Timing = Spike_TS(find(clust_idx==clust_no));
       % Spike_Timing = Spike_TS;
        
        figure
       size(Spike_Timing)
       hist(Spike_Timing,0:60:max(Spike_Timing));
        
        
        frame = zeros(size(Spike_Timing));
        for i=1:size(Spike_Timing,2)   % Future: a better way than this?
            frm = invoke(TTX, 'QryEpocAtV', 'fTrg', Spike_Timing(i), 0);
            frame(i) = frm;
            % The last param: 0 indicating the epoch's value returned;
            % 1 to get the epoch's start time or 2 to get the epoch's stop time.  
            % Pass 3 if you want the epocs filter status.
            if (frm>0) & (frm<=n_frames)
             n_spikes(frm)=n_spikes(frm)+1;
            end
        end


figure
bar(condensedata(n_spikes,15));
% 
% figure
% plot(1-cos(2*pi*(1:4500)/150));

%%% evaluate these before normalizing to rate
N=sum(n_spikes)
duration_wn(cell_n)=max(Spike_Timing);
sta_N(cell_n)=N;
fano(cell_n) = var(n_spikes)/mean(n_spikes)
    
    if N==0
        break
    end
n_spikes = n_spikes/(.03*n_reps);


    sta_responsiveness(cell_n) = 0;
    sta_contrastphase(cell_n) = 0;
    halfcontrast_wn(cell_n)=0;
    halfslope_wn(cell_n)=0;
    adaptation_index(cell_n)=0;


%%%plot histogram
    titlestr = sprintf('channel %d cluster %d',channel_no, clust_no);
%     figure
%     plot(n_spikes);
%     title(titlestr);
%     saveas(gcf,fullfile(pname,sprintf('hist%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');

for onoff=1:3 

    sta_length = 10;  %%% number of time points to calculate STA at
    sta_t = zeros(size(moviedata(:,:,1:sta_length)));
    N=0;
    sta2 = 0;    %%%spike triggered variance (not covariance!)

    tic
    for t = sta_length:n_frames
       if n_spikes(t)>0
            movie_snip =(moviedata(:,:,t-(sta_length-1):t));
            if onoff==1
                movie_snip(movie_snip==255)=128;
            elseif onoff==2
                movie_snip(movie_snip==0)=128;
           end
            sta_t =sta_t + n_spikes(t).*movie_snip;
            sta2 = sta2 + n_spikes(t).*movie_snip.^2;
            N = N + n_spikes(t);
       end
    end
    toc
    
 
    sta_t = sta_t(:,:,sta_length:-1:1)/N;
      sta2 = sta2(:,:,sta_length:-1:1)/N;
      stvar = (sta2 - sta_t.^2)/(128^2);
     


  
   if onoff==1
       color_range = [-12 12]
   else
       color_range = [-12 12]
   end
  
    %%% plot STA at each time point
    figure
    for t = 1:9
        subplot(3, 3, t);
        %imagesc(sta_t(:,:,t)'-movavg' ,color_range);
        imagesc(imfilter(sta_t(:,:,t)'-movavg',filt) ,color_range);
         sta_all(cell_n,:,:,t) = sta_t(:,:,t)-movavg;
        if t==1
            title_text = sprintf('%s',Tank_Name);
            text(0,-10,title_text,'FontSize',8);
        end
        if t==2
            title_text = sprintf('ch%d c%d',channel_no,clust_no);
            text(0,-10,title_text,'FontSize',8);  
        end
    end
    if printfig
        print(gcf);
    end

        if onoff==1
            saveas(gcf,fullfile(pname,sprintf('sta_off_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
        elseif onoff==2
              saveas(gcf,fullfile(pname,sprintf('sta_on_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
        else
                      saveas(gcf,fullfile(pname,sprintf('sta_sum_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
        end
                      %%% get subfield centered on RF and perform fourier analysis
    t=3;
    subplot(3,3,t);

    [x(1) x(2)] = (ginput(1));
    x=round(x)
    sta = get(gco,'CData');
    x= max(x,16);
    x= min(x,44);
   % subfield = sta(x(2)-23:x(2)+24,x(1)-23:x(1)+24);
       subfield = sta(x(2)-15:x(2)+16,x(1)-15:x(1)+16); 
       
        crop_sta = zeros(size(sta));
        crop_sta(x(2)-15:x(2)+16,x(1)-15:x(1)+16)=subfield;
 %   crop_sta = subfield;
%       x= max(x,16);
%     x= min(x,44);
%     subfield = sta(x(2)-15:x(2)+16,x(1)-15:x(1)+16);

    
    [sta_wpref(cell_n)  sta_thetapref(cell_n) sta_A1(cell_n) sta_thetawidth(cell_n) sta_baseline(cell_n) sta_null(cell_n)] ...
            = sta_analyze(crop_sta);
    

    %%% plot STV at each time point (this is some indication of complex cell response and RF location)
    figure
    for t = 1:9
        subplot(3, 3, t);
        imagesc(stvar(:,:,t)'-mov_var/(128^2));
    end
    title(titlestr);
        saveas(gcf,fullfile(pname,sprintf('stvar_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');

        figure
    for t = 1:9
        subplot(3, 3, t);
        imagesc(imfilter(sta_t(:,:,t)'-movavg',filt) ,color_range);
    end    
%         %%% time point of peak response
%     t_lag = 4;
%    
%     %%% fourier spectrum of STA
%     figure
%     imagesc(fftshift(abs(fft2(sta_t(:,:,t_lag)'-128))));
%      title(titlestr);
%     saveas(gcf,fullfile(pname,sprintf('fft_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
%      
%     h = zeros(size(n_spikes));
%     m=0;
%     tic
%     sta=(sta_t(xrange,yrange,t_lag)-movavg(xrange,yrange))/128;

end   %% on off
    
end    %%% cell


if use_afile
    save(afile,'sta_N','sta_wpref','sta_responsiveness','sta_contrastphase', 'halfslope_wn','response_wn','duration_wn','adaptation_index','sta_thetapref','sta_A1','sta_null','sta_thetawidth','sta_baseline','sta_all','transfer_function','fano','-append');
end
    
invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');