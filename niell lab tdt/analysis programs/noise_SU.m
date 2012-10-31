function noise_SU

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
for i =1:length(Block_Name);
    sprintf('%d : %s ',i,Block_Name{i})
end
block = input('which block to analyze ? ');
Block_Name = Block_Name{block}

load wn012alpha1_5hz_contrast10sec5min021407  %%% the usual
%load 5hz08cpdflat075sig10min
n_frames = 8999;
 %load 10hz10cpdflat075sig10min
 %n_frames =17999;
% load 012flat_nocmod_2to8Hz_05sig_10min
% n_frames = 17999;
% load sparse8deg5percent2x050108
% n_frames=17999
%load wn012cpd5hz30hz_alph1_off3_5min_021207
%load alpha_stepcontrast
%load sf_movie
clear amp width x0 baseline


contrast_modulated = 1;
correct_spectrum = 1;
frame_duration = 0.333;

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


calculate_stc =0;    %%% whether or not to calculate STC (it's slow!)
 xrange =25:45;
yrange = 30:50;   
dx = size(xrange,2);
dy=size(yrange,2);
moviedata = double(moviedata);
m = reshape(moviedata(xrange,yrange,1:n_frames),[dx*dy n_frames])/128;
if calculate_stc
    covmat_prior=0;


    tic
    covmat_prior = cov(double(m)');
    figure
    imagesc(covmat_prior);
    toc
end

if contrast_modulated
    for f = 1:n_frames;
         moviedata(:,:,f) = 128 + (moviedata(:,:,f)-128)/( 0.1 + 0.5 - 0.5*cos(2*pi*f/300));
    end
end

if correct_spectrum
    tic
    spectrum = zeros(size(moviedata),'single');
    for f= 1:n_frames
        spectrum(:,:,f) = fft2(single(moviedata(:,:,f))-128);
    end    
    mean_spectrum = mean(abs(spectrum(:,:,:)),3);
    mean_spectrum = mean_spectrum / max(max(mean_spectrum));
    toc

    inv_spectrum = (1./mean_spectrum).*(mean_spectrum>.05);

    tic
    for f = 1:n_frames
        moviedata(:,:,f) = ifft2(spectrum(:,:,f).*inv_spectrum)+128;
    end
    toc
end

m= m(:,1:1000);
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
  noisepname = pname;
  load(afile); 
  use_afile=1;
else
    dispay('need analysis file')
    return
end
%%% set time based on first and last timepoints in clustered data



for cell_n = 1:size(cells,1)
    channel_no = cells(cell_n,1)
    clust_no = cells(cell_n,2)
    channel_times =spikeT{cell_n} - (block-1)*10^5;
    times = channel_times(channel_times>0 & channel_times<10^5);
    
    sta0 = zeros(size(moviedata(:,:,1)));
    sta1 = sta0; sta2=sta0;
    frame_duration = median(diff(frameEpocs{block}(2,:)))
    
    n_spikes = zeros(n_frames,1);
    for f = 1:n_frames
        [Spike_Timing index numtrials] = getTrials(frameEpocs{block},times, f, frame_duration);
        n_spikes(f) = length(Spike_Timing);
        ntrials(f) = numtrials;
    end
        
n_reps = min(ntrials);
figure
bar(condenseData(n_spikes,15));
% 
% figure
% plot(1-cos(2*pi*(1:4500)/150));

%%% evaluate these before normalizing to rate
N=sum(n_spikes)
duration_wn(cell_n)=max(times);
sta_N(cell_n)=N;
fano(cell_n) = var(n_spikes)/mean(n_spikes)
    
    if N==0
        break
    end
n_spikes = n_spikes/(.03*n_reps);

cyc_frames = 10*30;
if contrast_modulated
%     for f= 1:n_frames;
%         frm = moviedata(:,:,f);
%         c(f) = double(std(frm(:)));
%     end
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
    figure
    loglog(abs(fftsignal));

    
    sta_responsiveness(cell_n) = 2*abs(fftsignal(31))/(abs(fftsignal(1)));  %%double to normalize FFT relative to DC
    sta_contrastphase(cell_n) = angle(fftsignal(31));
    phaseangle = angle(fftsignal(31))
    cycledata = zeros(cyc_frames,1);
    for i = 1:n_frames;
        cycledata(mod(i,cyc_frames)+1)=cycledata(mod(i,cyc_frames)+1)+n_spikes(i);
    end
    cycledata = cycledata./(nf/cyc_frames);
 
    %cycledata = conv(cycledata,ones(1,30))/10;
    cycledata = condenseData(cycledata,15);
    contrastdata = condenseData(contrastdata,15);
    
%      figure 
%     hold on  
%     plot(cycledata/max(cycledata),'b');
%     plot(double(contrastdata)/max(contrastdata),'g');

    %%% calculate contrast-response params
    contrastdata =double(contrastdata)/max(contrastdata);
    peakresp_wn(cell_n) = max(cycledata);
    norm_cycledata = (cycledata-min(cycledata(1:8)))/(max(cycledata(1:8))-min(cycledata(1:8)));
   %%% fit polynomial to rising side of contrast-response
   p = polyfit(double(contrastdata(1:8))/max(contrastdata),norm_cycledata(1:8),3);
%     figure
%      hold on
%     plot(double(contrastdata)/max(contrastdata),norm_cycledata); 
%     plot(0:0.1:contrastdata(8),polyval(p,double(0:0.1:contrastdata(8))),'go');
   
    %%% find contrast where response crosses 0.5, by finding polynomial root
    contrast_wn(:,cell_n) = contrastdata;
    response_wn(:,cell_n)= cycledata;
    p(4) = p(4)-0.5;
%     r = roots(p);
%     realroots = find(real(r)==r & r>0);
%     hf = min(r(realroots));
%     if isempty(hf)
%         hf = 1;
%     end
%     halfcontrast_wn(cell_n) = hf;
%     plot(halfcontrast_wn(cell_n),0.5,'r*');
%     halfslope_wn(cell_n) = polyval(polyder(p),halfcontrast_wn(cell_n));
    
    figure
    hold on   
   %%% plot average of up and down on contrast response
   meancontrast = 0.5*(contrastdata(1:10)+contrastdata(20:-1:11));
    meandata = 0.5*(cycledata(1:10)+cycledata(20:-1:11));
    plot(meancontrast/max(meancontrast),meandata/max(meandata),'g');
    title(sprintf('unit %d %d',channel_no, clust_no))
    xlabel('contrast');
    ylabel('response');
     plot(contrastdata/max(contrastdata),cycledata/max(meandata));
    saveas(gcf,fullfile(noisepname,sprintf('CRF_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
     %%% calculate adaptation index by difference between 
     %%%  up leg and down leg
   adaptation_index(cell_n) = sum(cycledata(1:10)-cycledata(20:-1:11))/sum(cycledata(1:10)+cycledata(20:-1:11));
     

else
    sta_responsiveness(cell_n) = 0;
    sta_contrastphase(cell_n) = 0;
    halfcontrast_wn(cell_n)=0;
    halfslope_wn(cell_n)=0;
    adaptation_index(cell_n)=0;
    
    
end

%%%plot histogram
    titlestr = sprintf('channel %d cluster %d',channel_no, clust_no);
%     figure
%     plot(n_spikes);
%     title(titlestr);
%     saveas(gcf,fullfile(pname,sprintf('hist%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');

 
    sta_length = 10;  %%% number of time points to calculate STA at
    sta_t = zeros(size(moviedata(:,:,1:sta_length)));
    N=0;
    sta2 = 0;    %%%spike triggered variance (not covariance!)

    tic
    for t = sta_length:n_frames
       if n_spikes(t)>0
            movie_snip =(moviedata(:,:,t-(sta_length-1):t));
            sta_t =sta_t + n_spikes(t).*movie_snip;
            sta2 = sta2 + n_spikes(t).*movie_snip.^2;
            N = N + n_spikes(t);
       end
    end
    toc
    
 
    sta_t = sta_t(:,:,sta_length:-1:1)/N;
      sta2 = sta2(:,:,sta_length:-1:1)/N;
      stvar = (sta2 - sta_t.^2)/(128^2);
     


  
   if correct_spectrum
       color_range = [-64 64]
   else
       color_range = [-64 64];
   end
  
    %%% plot STA at each time point
    figure
    for t = 1:9
        subplot(3, 3, t);
        imagesc(sta_t(:,:,t)'-movavg' ,color_range);
        axis square
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

  
    if contrast_modulated
        saveas(gcf,fullfile(noisepname,sprintf('sta_cmod_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
    else
         saveas(gcf,fullfile(noisepname,sprintf('sta_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
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
        saveas(gcf,fullfile(noisepname,sprintf('stvar_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');

        
        %%% time point of peak response
    t_lag = 4;
   
    %%% fourier spectrum of STA
%     figure
%     imagesc(fftshift(abs(fft2(sta_t(:,:,t_lag)'-128))));
%      title(titlestr);
%     saveas(gcf,fullfile(noisepname,sprintf('fft_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
     
    h = zeros(size(n_spikes));
    m=0;
    tic
    sta=(sta_t(xrange,yrange,t_lag)-movavg(xrange,yrange))/128;
    
    
    sta = sta/sqrt(sum(sum(sta.^2)));
    
    %%% calculate transfer function
    for t = (t_lag+1):n_frames
        if n_spikes(t)>=0
            m = (moviedata(xrange,yrange,t-t_lag+1)-movavg(xrange,yrange))/(128*sqrt(dx*dy));
            h(t) = sta(:)'*m(:);
        end
     end
    
%      figure
%      plot(h(t_lag+1:n_frames),n_spikes(t_lag+1:n_frames),'.');

     hist_int =0.05;
     h_round = round(h(t_lag+1:n_frames)/hist_int);
     n_sp = n_spikes(t_lag+1:n_frames);    
     h_min = min(h_round);
     h_range = h_min:max(h_round)-1
     n_mean=0; n_std=0; n_samp=0;
     for i = h_range;
         use = find(h_round ==i);
         n_mean(i-h_min+1) = mean(n_sp(use));
         n_std(i-h_min+1) = std(n_sp(use));
         n_samp(i-h_min+1) = size(use,1);
     end
     figure
     errorbar(h_range*hist_int,n_mean,n_std./sqrt(n_samp));
     hold on;
     plot(h_range*hist_int,n_samp/500,'g');

     transfer_function(cell_n,1:size(n_mean,2)) = n_mean;
     
    if calculate_stc
         covmat =0;
        N=0;
        figure
        imagesc(imresize(sta_t(xrange,yrange,t_lag)'-128,1,'bilinear'),[-32 32]);

        
        figure
        imagesc(stvar(xrange,yrange,t_lag)'-mov_var/(128^2),[-.1 .1]);
         %%% calculate spike triggered covariance matrix
         tic
         for t = (t_lag+1):n_frames
            if n_spikes(t)>0
                m = (moviedata(xrange,yrange,t-t_lag+1)-sta_t(xrange,yrange,t_lag))/128;
                covmat = covmat + n_spikes(t).*m(:)*m(:)';   
                N = N+n_spikes(t);
            end
         end
         toc

        covmat = covmat/N;
        figure
        imagesc(covmat);
        figure
        imagesc(covmat-covmat_prior,[-.05 .05])
        colorbar

        %%% find eigenvectors of the difference between spike-triggered and stimulus ensemble covariance 
        %%% this will find subspace spanned by RF as long as noise is gaussian (doesn't need to be white)
        [V D] = eig(covmat-covmat_prior);


        for i=1:dx*dy
            var_diag(i)=covmat(i,i)-covmat_prior(i,i);
        end
        stvar = reshape(var_diag,dx,dy);
        figure
        im = imresize(stvar,4,'bilinear');
        imagesc(im',[-.02 .04]);

        %%% sort and display eigenvalues
        clear lam V_sort
        for i=1:dx*dy; lam(i)=D(i,i); end
        figure
        [lam_sort I] = sort(lam);
        V_sort(:,:) = V(:,I);
         plot((lam_sort),'o');


        for i = dx*dy:-1:dx*dy-2
            im = reshape(V_sort(:,i),dx,dy);
            im = imresize(im,4,'bilinear');
            figure
            imagesc(im',[-0.2 0.2]);
            axis image
            if i==dx*dy
                    saveas(gcf,fullfile(noisepname,sprintf('stc_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
            end
        end

        %%% calculate transfer function for first eigenvector
        h = zeros(size(n_spikes));
        m=0;
        tic
         for t = (t_lag+1):n_frames
            if n_spikes(t)>=0
                m = (moviedata(xrange,yrange,t-t_lag+1)-sta_t(xrange,yrange,t_lag))/(128*sqrt(dx*dy));
                h(t) = V_sort(:,dx*dy)'*m(:);

            end
         end
         toc

         figure
         plot(h(t_lag+1:n_frames),n_spikes(t_lag+1:n_frames),'.');

         hist_int = .05;
         h_round = round(h(t_lag+1:n_frames)/hist_int);
         n_sp = n_spikes(t_lag+1:n_frames);    
         h_min = min(h_round)+1;
         h_range = h_min:(max(h_round)-1);
        n_mean=0; n_std=0; n_samp=0;
        for i = h_range;
             use = find(h_round ==i);
             n_mean(i-h_min+1) = mean(n_sp(use));
             n_std(i-h_min+1) = std(n_sp(use));
             n_samp(i-h_min+1) = size(use,1);
         end
         figure
         errorbar(h_range*hist_int,n_mean,n_std./sqrt(n_samp));
         hold on;
         plot(h_range*hist_int,n_samp/500,'g');
    end   %% if spike-triggered covariance 
    
end    %%% cell


if use_afile
    save(afile,'sta_N','sta_wpref','sta_responsiveness','sta_contrastphase', 'halfcontrast_wn','halfslope_wn','contrast_wn','response_wn','duration_wn','adaptation_index','sta_thetapref','sta_A1','sta_null','sta_thetawidth','sta_baseline','sta_all','transfer_function','fano','-append');
end
    
invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');