%function noise_analysis

% Matlab codes for reading from TTank for movie data
% Calculates RFs from spike triggered average, and spike triggered covariance
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03


%%% read in cluster data, then connect to the tank and read the block
clear all

SU = menu('recording type','multi-unit','single unit')-1
cells =1;
if SU
    [fname, pname] = uigetfile('*.mat','cluster data');
    load(fullfile(pname,fname));
    %     for i =1:length(Block_Name);
    %         sprintf('%d : %s ',i,Block_Name{i})
    %     end
    %block = input('which block to analyze ? ');
    block = listdlg('ListString',Block_Name,'SelectionMode','single')
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

[fname pname] = uigetfile('*.mat','movie file');
load(fullfile(pname,fname));

%load wn012alpha1_5hz_contrast10sec5min021407   %%% the usual

n_frames = length(moviedata)-1;

movietype = menu('Movie type','Contrast-modulated noise','Flashing sparse','Moving sparse');
cm_noise = 1;
fl_noise=2;
mv_noise=3;

if movietype==cm_noise
    contrast_modulated = input('contrast modulated? 0/1 : ');
    correct_spectrum = input('correct spectrum? 0/1 : ');
    pos_neg = input('separate on/off sta? 0/1 : ');
    compute_svd = input('compute svd ? 0/1 : ');
else
    contrast_modulated=0;
    correct_spectrum=0;
    pos_neg = input('separate on/off sta? 0/1 : ');
    compute_svd = input('compute svd ? 0/1 : ');
end


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

figure
imagesc(movavg-127,[-64 64])

if SU
    cell_range = 1:size(cells,1)
else
    cell_range=1:4:nchan
end
cell_range=5;
for cell_n = cell_range
    if SU
        channel_no = cells(cell_n,1)
        clust_no = cells(cell_n,2)
        channel_times =spikeT{cell_n} - (block-1)*10^5;
        times = channel_times(channel_times>0 & channel_times<10^5);
        frame_duration = median(diff(frameEpocs{block}(2,:)))
    else
        clust_no = [];
        channel_no = cell_n;
        frame_duration = median(diff(data.frameEpocs(2,:)))
        times=data.MUspikeT{cell_n};
    end
    if ~isempty(times)
        sta0 = zeros(size(moviedata(:,:,1)));
        sta1 = sta0; sta2=sta0;
        
        
        n_spikes = zeros(n_frames,1);
        for f = 1:n_frames
            if SU
                [Spike_Timing index numtrials] = getTrials(frameEpocs{block},times, f, frame_duration);
                eps = frameEpocs{block};
            else
                [Spike_Timing index numtrials] = getTrialsSU(data.frameEpocs,data.MUspikeT{cell_n}, f, frame_duration);
                eps = data.frameEpocs;
            end
            n_spikes(f) = length(Spike_Timing);
            ntrials(f) = numtrials;
        end
        
        n_reps = min(ntrials);
        

       
        
        %%% evaluate these before normalizing to rate
        N=sum(n_spikes)
        duration_wn(cell_n)=max(times);
        sta_N(cell_n)=N;
        fano(cell_n) = var(n_spikes)/mean(n_spikes);
        
        if N==0
            break
        end
        n_spikes = n_spikes./(frame_duration*ntrials');
       
       if movietype == cm_noise
            figure
        bar(condenseData(n_spikes,15));
        else
            figure
            hist(n_spikes)
        end
        
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
            if SU
                saveas(gcf,fullfile(noisepname,sprintf('CRF_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
            end     %%% calculate adaptation index by difference between
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
        %     saveas(gcf,fullfile (pname,sprintf('hist%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
        
        
        sta_length = 16;  %%% number of time points to calculate STA at
        sta_t = zeros(size(moviedata(:,:,1:sta_length)));
        N=0;
        sta2 = sta_t;    %%%spike triggered variance (not covariance!)
        sta_pos = sta_t;
        sta_neg = sta_t;
        tic
        for t = sta_length:n_frames
            if n_spikes(t)>0
                movie_snip =(moviedata(:,:,t-(sta_length-1):t));
                movie_snip_pos = (movie_snip-127)/127;
                movie_snip_neg =movie_snip_pos;
                movie_snip_pos(movie_snip_pos<0)=0;
                movie_snip_neg(movie_snip_pos>0)=0;
                sta_pos = sta_pos+n_spikes(t).*movie_snip_pos;
                sta_neg = sta_neg+n_spikes(t).*movie_snip_neg;
                sta_t =sta_t + n_spikes(t).*movie_snip;
                sta2 = sta2 + n_spikes(t).*((movie_snip-128).^2);
                N = N + n_spikes(t);
            end
        end
        toc
        
        
        sta_t = sta_t(:,:,sta_length:-1:1)/N;
        sta2 = sta2(:,:,sta_length:-1:1)/N;
        sta_pos = sta_pos(:,:,sta_length:-1:1)/N;
        sta_neg = sta_neg(:,:,sta_length:-1:1)/N;
        %stvar = (sta2 - sta_t.^2)/(128^2);
        stvar = sta2/(128^2);
        if movietype ==cm_noise
color_range = [-64 64];
        else
        color_range = [-32 32];
        end
        %%% plot STA at each time point
        figure
        for t = 1:16
            subplot(4, 4, t);
            imagesc(sta_t(:,:,t)'-movavg' ,color_range);
            axis equal
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
        
            if SU
            saveas(gcf,fullfile(noisepname,sprintf('sta_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
        end
        
        
        if movietype~=cm_noise
            sz_all = zeros(length(n_spikes),1);
            sz_rand = sz_all;
            sp_all = sz_all;
            th=sz_all;
            spd=sz_all;
            mov = sz_all;
            mov_rand= sz_all;
            [m ind] = max(abs(sta_t(:)-127));
            (m-128) / 128
            [x y t] = ind2sub(size(sta_t),ind)
            t_lag = t-1;
            
            movie_peak(cell_n,1)=x;
            movie_peak(cell_n,2)=y;
            movie_lag(cell_n) = t_lag*frame_duration;
            h=0;
            h_rand=0;
            
            mh = 0;
            mh_rand=0;
            mhistbins = [0 127 255];
             histbins = [0 1 2 4 8 16];
%             
%             if movietype == fl_noise
%                 histbins = [0 1 2 4 8 16];
%             elseif movietype == mv_noise
%                 histbins = [0 1 2 3 4 5];
%             end
            w=0;
            x = max(x,w+1);
            y= max(y,w+1);
            x= min(x,size(moviedata,1)-w);
            y= min(y,size(moviedata,2)-w);
            spd_hist=0; spd_hist_rand=0;
            th_hist=0; th_hist_rand=0;
            spdbins = [0 10 20 40 80 160];
            thbins = [0:pi/4:2*pi];
            
            
            for f = t_lag+1:length(n_spikes);
                mov(f) = double(moviedata(x,y,f-t_lag));
                mov_rand(f) =double(moviedata(round(rand)*127+1,round(rand*127)+1,f-t_lag));
                sz_all(f) = double(sz_mov(x,y,f-t_lag));
                sz_rand(f) = double(sz_mov(round(rand)*127+1,round(rand*127)+1,f-t_lag));
                
                sp_all(f) = n_spikes(f);
                if movietype == mv_noise
                    spd(f) = double(sp_mov(x,y,f-t_lag));
                    th(f) =double(th_mov(x,y,f-t_lag))*2*pi/255;
                end
                if sp_all(f)>0
                    if movietype ==mv_noise
                        spd_hist = spd_hist+hist(double(sp_mov(x-w:x+w,y-w:y+w,f-t_lag)),spdbins)*sp_all(f);
                        th_hist = th_hist+hist(double(th_mov(x-w:x+w,y-w:y+w,f-t_lag))*2*pi/255,thbins)*sp_all(f);
                        spd_hist_rand = spd_hist_rand+sp_all(f)*hist(double(sp_mov(x-w:x+w,y-w:y+w,ceil(rand*length(n_spikes)))),spdbins);
                        th_hist_rand = th_hist_rand+sp_all(f)*hist(double(th_mov(x-w:x+w,y-w:y+w,ceil(rand*length(n_spikes))))*2*pi/255,thbins);
                    end
                    
                    h = h+hist(double(sz_mov(x-w:x+w,y-w:y+w,f-t_lag)),histbins)*sp_all(f);
                    h_rand = h_rand+hist(double(sz_mov(x-w:x+w,y-w:y+w,ceil(rand*length(n_spikes)))),histbins)*sp_all(f);
                    
                    mh = mh+hist(mov(f),mhistbins)*sp_all(f);
                    mh_rand = mh_rand+hist(mov_rand(f),mhistbins)*sp_all(f);
                end
            end
            
            %         figure
            %         polar(th,spd,'o');
            %         hold on
            %         polar(th(sp_all>0),spd(sp_all>0),'go');
            
            if movietype == mv_noise
                sp_norm = spd_hist/sum(sp_all);
                sp_norm_rand = spd_hist_rand/sum(sp_all);
                
                figure
                plot(sp_norm);
                hold on
                plot(sp_norm_rand,'r');
                xlabel('speed')
                
                th_norm = th_hist/sum(sp_all);
                th_norm_rand = th_hist_rand/sum(sp_all);
                figure
                plot(th_norm);
                hold on
                plot(th_norm_rand,'r');
                xlabel('theta')
                
                figure
                plot((sp_norm(2:6)-sp_norm_rand(2:6))./sp_norm_rand(2:6))
                xlabel('speed')
                spd_tuning(cell_n,:) = (sp_norm(2:6)-sp_norm_rand(2:6))./sp_norm_rand(2:6);
                %         figure
                %         plot(th_hist);
                %
            end
            
            figure
            plot(h/nansum(sp_all));
            hold on
            plot(h_rand/nansum(sp_all),'r');
            xlabel('size')
            
            sz_tune_rand = h_rand/nansum(sp_all);
            sz_tune = h/nansum(sp_all);
            baseline = (1-sz_tune_rand(1))/(length(sz_tune)-1)
            
            figure
            plot((sz_tune(2:6)-baseline)/baseline)
            xlabel('size')
            
            size_tuning(cell_n,:) =(sz_tune(2:6)-baseline)/baseline;
            
            %%%%% histrograms relative to onset/offset
            rastfig = figure;
            histfig=figure;
            
            for rep = 1:2
                n=0;
                ontime=0;
                tseries = squeeze(moviedata(x,y,:));
                if rep ==1
                    onset = (tseries(1:end-1)==127 & diff(tseries)>0);
                    c = 'r';
                elseif rep ==2
                    onset = (tseries(1:end-1)==127 & diff(tseries)<0);
                    c= 'b';
                elseif rep ==3
                    onset = (tseries(1:end-1)==255 & diff(tseries)<0);
                    c= 'g';
                elseif rep ==4
                    onset = (tseries(1:end-1)==0 & diff(tseries)>1);
                    c= 'c';
                end
                
                for i=1:length(eps);
                    if eps(1,i)<length(onset) & onset(eps(1,i))
                        n= n+1;
                        ontime(n) = eps(2,i);
                    end
                end
                
                
                if movietype == fl_noise
                    dt=0.05
                    histbins = 0:dt:0.75;
                elseif movietype==mv_noise
                    dt = 0.1
                    histbins = -1:dt:2;
                end
                
                h=0;
                figure(rastfig)
                hold on
                for i = 1:n
                    t = times-ontime(i);
                    t = t(t>histbins(1) & t<histbins(end));
                    plot(t,ones(length(t),1)*i,[c '*'])
                    if ~isempty(t)
                        h =h+ histc(t, histbins);
                    end
                end
                
               if length(h)>1
                   figure(histfig)
                hold on
                plot(histbins(1:end-1)+dt/2,h(1:end-1)/(n*dt),c)
                onset_hist(cell_n,rep,:) = h(1:end-1)/(n*dt);
               end
            end
            figure(histfig);
            legend('gray->on','gray->off')
            if movietype ==fl_noise
                plot([0.25 0.5], [0.5 0.5],'g','LineWidth',8)
            end
            
            %         figure
            %         plot(mh/nansum(sp_all));
            %         hold on
            %         plot(mh_rand/nansum(sp_all),'r');
            %
            
            
        end
        nx = size(sta_t,1);
        ny = size(sta_t,2);
        nt = size(sta_t,3);
        
        
        
        if compute_svd
            sta_col =reshape(sta_t,nx*ny,nt);
            [u s v] = svd(sta_col-mean(sta_t(:)));
            
            figure
            
            for i = 1:3
                subplot(1,3,i)
                range = max(abs(min(u(:,i))),abs(max(u(:,i))));
                imagesc(reshape(u(:,i),nx,ny)'*sign(v(1,i)),1.2*[-range range]);
                v(:,i) = v(:,i)*sign(v(1,i));
                axis square;
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
            end
            figure
            plot(v(:,1:3));
            figure
            plot(diag(s));
        end
        
        if pos_neg
            color_range = [0 1];
            figure
            for t = 1:9
                subplot(3, 3, t);
                imagesc(sta_pos(:,:,t)' ,color_range);
                axis equal
                
                if t==1
                    title_text = sprintf('%s',Tank_Name);
                    text(0,-10,title_text,'FontSize',8);
                end
                if t==2
                    title_text = sprintf('ch%d c%d',channel_no,clust_no);
                    text(0,-10,title_text,'FontSize',8);
                end
            end
            
            figure
            color_range = [-1 0];
            for t = 1:9
                subplot(3, 3, t);
                imagesc(sta_neg(:,:,t)' ,color_range);
                axis equal
                
                if t==1
                    title_text = sprintf('%s',Tank_Name);
                    text(0,-10,title_text,'FontSize',8);
                end
                if t==2
                    title_text = sprintf('ch%d c%d',channel_no,clust_no);
                    text(0,-10,title_text,'FontSize',8);
                end
            end
        end
        
        %     if SU
        %         if contrast_modulated
        %             saveas(gcf,fullfile(noisepname,sprintf('sta_cmod_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
        %         else
        %             saveas(gcf,fullfile(noisepname,sprintf('sta_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
        %         end
        %     end
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
        %
        %
        %     %    [sta_wpref(cell_n)  sta_thetapref(cell_n) sta_A1(cell_n) sta_thetawidth(cell_n) sta_baseline(cell_n) sta_null(cell_n)] ...
        %     %        = sta_analyze(crop_sta);
        %
        %
        %
        %     %%% plot STV at each time point (this is some indication of complex cell response and RF location)
        %     figure
        %     for t = 1:9
        %         subplot(3, 3, t);
        %         imagesc(stvar(:,:,t)'-mov_var/(128^2),[-0.1 0.1]);
        %     end
        %     title(titlestr);
        
        if SU
            saveas(gcf,fullfile(noisepname,sprintf('stvar_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
        end
        
        %%% time point of peak response
        %t_lag = 4;
        
        %% fourier spectrum of STA
        if movietype==cm_noise
            [m ind] = max(abs(sta_t(:)-127));
            
            [x y t] = ind2sub(size(sta_t),ind)
            figure
            imagesc(fftshift(abs(fft2(sta_t(:,:,t)'-128))));
             title(titlestr);
            saveas(gcf,fullfile(noisepname,sprintf('fft_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
        end
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
        h_range = h_min:max(h_round)-1;
        n_mean=0; n_std=0; n_samp=0;
        for i = h_range;
            use = find(h_round ==i);
            n_mean(i-h_min+1) = mean(n_sp(use));
            n_std(i-h_min+1) = std(n_sp(use));
            n_samp(i-h_min+1) = size(use,1);
        end
        if movietype == cm_noise
            figure
        errorbar(h_range*hist_int,n_mean,n_std./sqrt(n_samp));
        hold on;
        plot(h_range*hist_int,n_samp/500,'g');
        end
        
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
    end
end    %%% cell

% 
% if SU
%     save(afile,'sta_N','sta_wpref','sta_responsiveness','sta_contrastphase','halfslope_wn','contrast_wn','response_wn','duration_wn','adaptation_index','sta_thetapref','sta_A1','sta_null','sta_thetawidth','sta_baseline','sta_all','transfer_function','fano','-append');
% end
