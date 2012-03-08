%function noise_cluster

% Matlab codes for reading from TTank for movie data
% Calculates RFs from spike triggered average, and spike triggered covariance
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03


%%% read in cluster data, then connect to the tank and read the block
% clear all
% [fname, pname] = uigetfile('*.mat','cluster data');
% oldpname = pname;
% load(fullfile(pname,fname));
% pname = oldpname;
% 
% % generally analyze cluster data that is bars/drift/movie (3) or movie
% % alone (1)
% if nblocks>1
%     block=2;
%     Block_Name = char(Block_Name(block)) % to initialize
% else
%     block=1;
%     Block_Name = char(Block_Name(block)) % to initialize - darcy_90616
% end

%%% read in cluster data, then connect to the tank and read the block
cells=0;
[fname, pname] = uigetfile('*.mat','cluster data');
oldpname = pname;
load(fullfile(pname,fname));   %%% need to copy pname, or it can get written over in load
pname = oldpname;
for i =1:length(Block_Name);
    sprintf('%d : %s ',i,Block_Name{i})
end
block = input('which block to analyze ? ');
Block_Name = Block_Name{block}
TTX = openTTX(Tank_Name,Block_Name); % to initialize



load 05deg4hzmovie %%% movie frames
n_frames = 1799;

clear amp width x0 baseline
Tank_Name
TTX = openTTX(Tank_Name,Block_Name); % to initialize

contrast_modulated = 1;
correct_spectrum = 0;

Event_Name_Snip='Snip'
Event_Name_Wave='PDec'
Sample_Interval=0.04096; % 24414.0625Hz
Sample_Number_Snip=64;
Dec_Factor=32;
plot_duration=36; %in second
hist_int = 0.1;
hist_range=[0:hist_int:plot_duration];
axis_range=[0 plot_duration 0 20];
max_events=100000;

xrange = 25:45;
yrange = 30:50;
dx = size(xrange,2);
dy = size(yrange,2);
moviedata = double(moviedata);
m = reshape(moviedata(xrange,yrange,1:n_frames),[dx*dy n_frames])/128;

if contrast_modulated
    for f = 1:n_frames;
        moviedata(:,:,f) = 128 + (moviedata(:,:,f)-128)/( 0.1 + 0.5 - 0.5*cos(2*pi*f/300));
    end
end

m = m(:,1:1000);
mov_var = var(squeeze(moviedata(50,50,:)))
movavg = mean(moviedata,3);

tetrode_linear=0;
if tetrode_linear
    ch_map = [14 8 10 4 13 7 9 3 11 1 15 5 12 2 16 6];
else
    ch_map = 1:16;
end
[fname, pname] = uigetfile('*.mat','analysis data');
if fname~=0
    afile = fullfile(pname,fname);
    load(afile);
    use_afile=1;
else
    use_afile=0;
    %%% select which units to analyze (channel, cluster number)
    cells = [7 1; 10 1; 12 1; 12 2; 14 1; 15 2 ; 15 3; 16 2]
    cells = [1 4; 5 1; 9 1; 9 6; 9 2; 13 2; 13 3; 13 5];
    cells = [1 5; 5 4; 5 6; 9 2; 9 4; 9 6; 13 2; 13 4; 13 6; 13 7 ]
    cells = [1 4; 5 1; 9 1; 9 2; 9 6; 13 1; 13 4; 13 6; 13 7];
    cells = [1 4; 5 4; 5 7; 9 2; 9 7; 13 3; 13 5; 13 7]  %%% 011606_allstim1
    cells = [1 4; 5 4; 5 6; 9 4; 9 6; 13 1; 13 4; ]
    cells = [1 2; 1 4; 1 7; 9 1; 9 2; 9 4; 9 6; 13 1; 13 5; ]  %%% 011606_allstim2
    cells = [1 1;  5 4; 9 4; 9 5; 13 1]
    cells = [1 1; 5 1; 5 2; 5 7; 9 4; 9 5;  9 3; 13 4]; %%%011707_allstim1
    cells = [5 4; 5 6];
    cells = [1 7; 5 7; 9 2; 13 7; 13 1] %% allstim 011806
    cells = [1 1; 1 2; 5 1; 5 2; 5 4; 5 6; 5 7; 9 4; 9 5;  9 3; 9 7; 13 4; 13 3; 13 7]; %%%011707_allstim1
end
%%% set time based on first and last timepoints in clustered data
invoke(TTX,'CreateEpocIndexing')

event_times_all = event_times_all-(block-1)*10^5;
times = event_times_all(event_times_all>0 & event_times_all<10^5);
time1 = min(times);
time2= max(times);

MyEpocs = invoke(TTX, 'GetEpocsV', 'fTrg',time1, time2, round(max_time*30));

cells

for cell_n = 1:size(cells,1)
    cell_n
    channel_no = cells(cell_n,1)
    clust_no = cells(cell_n,2)

    title_text = sprintf('%s ch%d c%d',Tank_Name, channel_no, clust_no);

    
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

%     figure
%     size(Spike_Timing)
%     hist(Spike_Timing,0:60:time2);


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
    % figure
    % plot(frame, Spike_Timing, '.k', 'MarkerSize',6);

    % figure
    % plot(n_spikes);

    hist_fig=figure('Name', sprintf('Spike histogram ch:%d c%d', channel_no, clust_no));
    bar(condensedata(n_spikes,15));
    title(title_text,'FontSize',8);

    %
%     figure
%     plot(1-cos(2*pi*(1:4500)/150));

    %%% evaluate these before normalizing to rate
    N=sum(n_spikes)
    duration_wn(cell_n)=time2;
    sta_N(cell_n)=N;
    fano(cell_n) = var(n_spikes)/mean(n_spikes)

    % if N==0         %  active in Dan's version, but causes analysis to
                       % halt if there is a cell of low spiking; therefore
                       % I inactivated this step
     %   break
    % end

%     normalize # of spikes
    n_spikes = n_spikes/(.03*n_reps); % .03 s = frame length

    cyc_frames = 10*30;
    if contrast_modulated
        f = 1:n_frames;
        c = double(0.5- 0.5*cos(2*pi*f/cyc_frames));  %%% define contrast

        nf = zeros(300,1);
        contrastdata = zeros(cyc_frames,1);
        for i = 1:n_frames;
            contrastdata(mod(i,cyc_frames)+1)=contrastdata(mod(i,cyc_frames)+1)+double(c(i));
            nf(mod(i,cyc_frames)+1) = nf(mod(i,cyc_frames)+1)+1;
        end
        contrastdata = contrastdata./(nf/cyc_frames);
        fftsignal = (fft(n_spikes));
        fftsignal = fftsignal(1:round(n_frames/2));

        % normalize fftsignal to total # of spikes
        fftsignal=fftsignal/length(n_spikes);
        fftsignal(2:length(fftsignal))=fftsignal(2:length(fftsignal))*2; % normalize for +/- component of fft components > 1
        
        fft_fig=figure('Name', sprintf('fft ch:%d c%d', channel_no, clust_no));
        loglog(abs(fftsignal));
        title(title_text,'FontSize',8);

     
        
        %         figure
%         plot(abs(fftsignal));

        sta_responsiveness(cell_n) = 2*abs(fftsignal(7))/(abs(fftsignal(1)));  %%double to normalize FFT relative to DC
        sta_contrastphase(cell_n) = angle(fftsignal(7));
  
        r_f0(cell_n) = abs(fftsignal(1));
        r_f1(cell_n) = 2*abs(fftsignal(7));
        r_f1_phase(cell_n)=angle(fftsignal(7));
        r_f1dividef0(cell_n) = 2*abs(fftsignal(7))/(abs(fftsignal(1)));  %%double to normalize FFT relative to DC
        
%         cycledata = zeros(cyc_frames,1);
%         for i = 1:n_frames;
%             cycledata(mod(i,cyc_frames)+1)=cycledata(mod(i,cyc_frames)+1)+n_spikes(i);
%         end
%         cycledata = cycledata./(nf/cyc_frames);        
%         %cycledata = conv(cycledata,ones(1,30))/10;
%         cycledata = condenseData(cycledata,15);
%         contrastdata = condenseData(contrastdata,15);
% 
%         figure
%         hold on
%         plot(cycledata/max(cycledata),'b');
%         plot(double(contrastdata)/max(contrastdata),'g');
% 
%         %%% calculate contrast-response params
%         contrastdata =double(contrastdata)/max(contrastdata);
%         peakresp_wn(cell_n) = max(cycledata);
%         norm_cycledata = (cycledata-min(cycledata(1:8)))/(max(cycledata(1:8))-min(cycledata(1:8)));
%         %%% fit polynomial to rising side of contrast-response
%         p = polyfit(double(contrastdata(1:8))/max(contrastdata),norm_cycledata(1:8),3);
%         figure
%         hold on
%         plot(double(contrastdata)/max(contrastdata),norm_cycledata);
%         plot(0:0.1:contrastdata(8),polyval(p,double(0:0.1:contrastdata(8))),'go');

%         %%% find contrast where response crosses 0.5, by finding polynomial root
%         contrast_wn(:,cell_n) = contrastdata;
%         response_wn(:,cell_n)= cycledata;
%         p(4) = p(4)-0.5;
%         r = roots(p);
%         realroots = find(real(r)==r & r>0);
%         hf = min(r(realroots));
%         if isempty(hf)
%             hf = 1;
%         end
%         halfcontrast_wn(cell_n) = hf;
%         plot(halfcontrast_wn(cell_n),0.5,'r*');
%         halfslope_wn(cell_n) = polyval(polyder(p),halfcontrast_wn(cell_n));
% 
%         figure
%         hold on
%         %%% plot average of up and down on contrast response
%         meancontrast = 0.5*(contrastdata(1:10)+contrastdata(20:-1:11));
%         meandata = 0.5*(cycledata(1:10)+cycledata(20:-1:11));
%         plot(meancontrast/max(meancontrast),meandata/max(meandata),'g');
%         title(sprintf('unit %d %d',channel_no, clust_no))
%         xlabel('contrast');
%         ylabel('response');
%         plot(contrastdata/max(contrastdata),cycledata/max(meandata));
% 
%         %%% calculate adaptation index by difference between
%         %%%  up leg and down leg
%         adaptation_index(cell_n) = sum(cycledata(1:10)-cycledata(20:-1:11))/sum(cycledata(1:10)+cycledata(20:-1:11));
% 
    else
        sta_responsiveness(cell_n) = 0;
        sta_contrastphase(cell_n) = 0;
        halfcontrast_wn(cell_n)=0;
        halfslope_wn(cell_n)=0;
        adaptation_index(cell_n)=0;
    end

%     %%%plot histogram
%     titlestr = sprintf('channel %d cluster %d',channel_no, clust_no);
%     %     figure
%     %     plot(n_spikes);
%     %     title(titlestr);
%     %     saveas(gcf,fullfile(pname,sprintf('hist%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
% 
% 
%     sta_length = 10;  %%% number of time points to calculate STA at
%     sta_t = zeros(size(moviedata(:,:,1:sta_length)));
%     N=0;
%     sta2 = 0;    %%%spike triggered variance (not covariance!)
% 
%     tic
%     for t = sta_length:n_frames
%         if n_spikes(t)>0
%             movie_snip =(moviedata(:,:,t-(sta_length-1):t));
%             sta_t =sta_t + n_spikes(t).*movie_snip;
%             sta2 = sta2 + n_spikes(t).*movie_snip.^2;
%             N = N + n_spikes(t);
%         end
%     end
%     toc
% 
%     sta_t = sta_t(:,:,sta_length:-1:1)/N;
%     sta2 = sta2(:,:,sta_length:-1:1)/N;
%     stvar = (sta2 - sta_t.^2)/(128^2);
% 
%     if correct_spectrum
%         color_range = [-64 64]
%     else
%         color_range = [-64 64];
%     end
% 
%     %%% plot STA at each time point
%     figure
%     for t = 1:9
%         subplot(3, 3, t);
%         imagesc(sta_t(:,:,t)'-movavg' ,color_range);
%         sta_all(cell_n,:,:,t) = sta_t(:,:,t)-movavg;
%         if t==1
%             title_text = sprintf('%s',Tank_Name);
%             text(0,-10,title_text,'FontSize',8);
%         end
%         if t==2
%             title_text = sprintf('ch%d c%d',channel_no,clust_no);
%             text(0,-10,title_text,'FontSize',8);
%         end
%     end
% 
%     if contrast_modulated
%         saveas(gcf,fullfile(pname,sprintf('sta_cmod_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
%     else
%         saveas(gcf,fullfile(pname,sprintf('sta_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
%     end
% 
%     %%% get subfield centered on RF and perform fourier analysis
%     t=3;
%     subplot(3,3,t);
% 
%     [x(1) x(2)] = (ginput(1));
%     x=round(x)
%     sta = get(gco,'CData');
%     x= max(x,16);
%     x= min(x,44);
%     % subfield = sta(x(2)-23:x(2)+24,x(1)-23:x(1)+24);
%     subfield = sta(x(2)-15:x(2)+16,x(1)-15:x(1)+16);
% 
%     crop_sta = zeros(size(sta));
%     crop_sta(x(2)-15:x(2)+16,x(1)-15:x(1)+16)=subfield;
%     %   crop_sta = subfield;
%     %       x= max(x,16);
%     %     x= min(x,44);
%     %     subfield = sta(x(2)-15:x(2)+16,x(1)-15:x(1)+16);
%     figure
%     imagesc(subfield,[-64 64]);
% 
%     [sta_wpref(cell_n)  sta_thetapref(cell_n) sta_A1(cell_n) sta_thetawidth(cell_n) sta_baseline(cell_n) sta_null(cell_n)] ...
%         = sta_analyze(crop_sta);
% 
%     figure
%     imagesc(sta_t(:,:,3)'-movavg' ,1.25*color_range);
%     axis equal
% 
%     %%% plot STV at each time point (this is some indication of complex cell response and RF location)
%     figure
%     for t = 1:9
%         subplot(3, 3, t);
%         imagesc(stvar(:,:,t)'-mov_var/(128^2));
%     end
%     title(titlestr);
%     saveas(gcf,fullfile(pname,sprintf('stvar_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
% 
%     figure
%     imagesc(stvar(:,:,2)'-mov_var/(128^2));
% 
%     %%% time point of peak response
%     t_lag = 4;
% 
%     %%% fourier spectrum of STA
%     figure
%     imagesc(fftshift(abs(fft2(sta_t(:,:,t_lag)'-128))));
%     title(titlestr);
%     saveas(gcf,fullfile(pname,sprintf('fft_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
% 
%     h = zeros(size(n_spikes));
%     m=0;
%     tic
%     sta=(sta_t(xrange,yrange,t_lag)-movavg(xrange,yrange))/128;
% 
% 
%     sta = sta/sqrt(sum(sum(sta.^2)));
% 
%     %%% calculate transfer function
%     for t = (t_lag+1):n_frames
%         if n_spikes(t)>=0
%             m = (moviedata(xrange,yrange,t-t_lag+1)-movavg(xrange,yrange))/(128*sqrt(dx*dy));
%             h(t) = sta(:)'*m(:);
%         end
%     end
% 
%     figure
%     plot(h(t_lag+1:n_frames),n_spikes(t_lag+1:n_frames),'.');
% 
%     hist_int =0.05;
%     h_round = round(h(t_lag+1:n_frames)/hist_int);
%     n_sp = n_spikes(t_lag+1:n_frames);
%     h_min = min(h_round);
%     h_range = h_min:max(h_round)-1
%     n_mean=0; n_std=0; n_samp=0;
%     for i = h_range;
%         use = find(h_round ==i);
%         n_mean(i-h_min+1) = mean(n_sp(use));
%         n_std(i-h_min+1) = std(n_sp(use));
%         n_samp(i-h_min+1) = size(use,1);
%     end
%     figure
%     errorbar(h_range*hist_int,n_mean,n_std./sqrt(n_samp));
%     hold on;
%     plot(h_range*hist_int,n_samp/500,'g');
% 
%     transfer_function(cell_n,1:size(n_mean,2)) = n_mean;

saveas(hist_fig,fullfile(pname,sprintf('moviehist_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
saveas(fft_fig,fullfile(pname,sprintf('moviefft_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');



end    %%% cell


% if use_afile
%     save(afile,'sta_N','sta_wpref','sta_responsiveness','sta_contrastphase', 'halfcontrast_wn','halfslope_wn','contrast_wn','response_wn','duration_wn','adaptation_index','sta_thetapref','sta_A1','sta_null','sta_thetawidth','sta_baseline','sta_all','transfer_function','fano','-append');
% end


if use_afile
    save(afile,'r_f0','r_f1','r_f1_phase','r_f1dividef0','-append');
end



invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');