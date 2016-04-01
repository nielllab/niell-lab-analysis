function noise_analysis(clustfile,afile,pdfFile,movieFile,Block_Name,blocknum,movietype, stim_eye,sess)
% Matlab codes for reading from TTank for movie data
% Calculates RFs from spike triggered average, and spike triggered covariance
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03


%%% read in cluster data, then connect to the tank and read the block
close all
dbstop if error

if ~exist('Block_Name','var');
    SU = menu('recording type','multi-unit','single unit')-1;
    useArgin=0;
else
    SU=1;
    useArgin=1;
end


cells =1;
if SU
    if ~useArgin
        [fname, pname] = uigetfile('*.mat','cluster data');
        clustfile=fullfile(pname,fname);
    end
    load(clustfile);
    if ~useArgin
        blocknum = listdlg('ListString',Block_Name,'SelectionMode','single');
        
        
        [afname, apname] = uigetfile('*.mat','analysis data');
        noisepname = apname;
        afile = fullfile(apname,afname);
    end
    load(afile);
    afile
    [noisepname noisefname] = fileparts(afile);
    Block_Name = Block_Name{blocknum}
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

Block_Name

if ~useArgin
    [fname pname] = uigetfile('*.mat','movie file','C:\Users\lab\Desktop\movie files\cortex\');
    movieFile = fullfile(pname,fname);
end
load(movieFile);

%load wn012alpha1_5hz_contrast10sec5min021407   %%% the usual

%n_frames = length(moviedata)-1;
n_frames = max(frameEpocs{blocknum}(1,:))
moviedata=moviedata(:,:,1:n_frames);
n_frames = n_frames-1; %%% last frame is sometimes bad

if ~useArgin
    movietype = menu('Movie type','Contrast-modulated noise','Flashing sparse','Moving sparse');
end
cm_noise = 1;
fl_noise=2;
mv_noise=3;

if movietype==cm_noise
    prompt = {'contrast modulated','frame rate','correct spectrum','crop to screen size','contra(1) or ipsi(2) eye'};
    num_lines = 1;
    def = {'1','30','1','1','1'};
    if ~useArgin
        answer = inputdlg(prompt,'wn parameters',num_lines,def);
    else
        answer=def;
    end
    contrast_modulated = str2num(answer{1})
    framerate = str2num(answer{2})
    correct_spectrum = str2num(answer{3})
    crop_mov =  str2num(answer{4})
    if ~useArgin
        stim_eye =  str2num(answer{5})
    end
    
    pos_neg=0;
    compute_svd=1;
    contrast_period=10;
else
    contrast_modulated=0;
    correct_spectrum=0;
    pos_neg=0;
    compute_svd=0;
    crop_mov=0;
    % pos_neg = input('separate on/off sta? 0/1 : ');
    %compute_svd = input('compute svd ? 0/1 : ');
end


if useArgin
    psfilename = [pdfFile(1:end-4) Block_Name '.ps'];
else
    [fname pname] =uiputfile('*.ps'); psfilename=fullfile(pname,fname);
end
if exist(psfilename,'file')==2;delete(psfilename);end



calculate_stc =0;    %%% whether or not to calculate STC (it's slow!)
xrange =25:45;
yrange = 30:50;
dx = size(xrange,2);
dy=size(yrange,2);
if movietype==cm_noise
    moviedata = double(moviedata);
end
%m = reshape(moviedata(xrange,yrange,1:n_frames),[dx*dy n_frames])/128;
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
        moviedata(:,:,f) = 128 + (moviedata(:,:,f)-128)/( 0.1 + 0.5 - 0.5*cos(2*pi*f/(framerate*contrast_period)));
    end
end

%%%
spectrum_thresh=0.2; %%%to avoid bringing up high freqs
%spectrum_thresh=0.05 %%% to correct all freqs, but not rounding error

if correct_spectrum
    tic
    spectrum = zeros(size(moviedata),'single');
    for f= 1:n_frames
        spectrum(:,:,f) = fft2(single(moviedata(:,:,f))-128);
    end
    mean_spectrum = mean(abs(spectrum(:,:,:)),3);
    mean_spectrum = mean_spectrum / max(max(mean_spectrum));
    toc
    figure
    imagesc(fftshift(mean_spectrum))
    inv_spectrum = (1./mean_spectrum).*(mean_spectrum>spectrum_thresh);
    
    tic
    for f = 1:n_frames
        moviedata(:,:,f) = ifft2(spectrum(:,:,f).*inv_spectrum)+128;
    end
    toc
end

if crop_mov
    moviedata=moviedata(:,27:102,:);
end


% figure
% imagesc(movavg-127,[-64 64])

parpool
tic

if movietype==mv_noise
    th_mov=single(th_mov);
    
    
    th_mov=th_mov*2*pi/255;
    th_mov(th_mov>(15*pi/8)) = th_mov(th_mov>(15*pi/8))-2*pi;
    th_mov= uint8(round((th_mov)/(pi/4)+1));
    th_mov(sz_mov==0)=255;
end
toc
%moviedata=single(moviedata);

sta_length = 16;contrast_period=10;
display('reshaping movie')
tic

mov_squeeze = reshape(moviedata(:,:,1:end-sta_length),size(moviedata,1)*size(moviedata,2),size(moviedata,3)-sta_length);
if movietype~=cm_noise
    mov_squeeze=uint8(mov_squeeze);
end
toc

display('computing avgs')
tic

if movietype==fl_noise | movietype==mv_noise
    mov_neg = (mov_squeeze==0);
    mov_zero = (mov_squeeze==127);
    mov_pos= (mov_squeeze==255);
    s = ones(size(mov_squeeze,2),sta_length,'single');
    mov_neg_avg =  reshape(single(mov_neg*s),size(moviedata,1),size(moviedata,2),sta_length)/length(s);
    mov_zero_avg =  reshape(single(mov_zero*s),size(moviedata,1),size(moviedata,2),sta_length)/length(s);
    mov_pos_avg =  reshape(single(mov_pos*s),size(moviedata,1),size(moviedata,2),sta_length)/length(s);
    
end
toc

clear wn

if SU
    cell_range = 1:size(cells,1)
else
    cell_range=1:4:nchan
end
%cell_range=5;

for cell_n = cell_range
    tic
    display('reading data')
    if SU
        channel_no = cells(cell_n,1)
        clust_no = cells(cell_n,2)
        channel_times =spikeT{cell_n} - (blocknum-1)*10^5;
        times = channel_times(channel_times>0 & channel_times<10^5);
        frame_duration = median(diff(frameEpocs{blocknum}(2,:)))
    else
        clust_no = [];
        channel_no = cell_n;
        frame_duration = median(diff(data.frameEpocs(2,:)))
        times=data.MUspikeT{cell_n};
    end
    toc
    
    if isempty(times)
        display('no spikes ... skipping')
    else
        n_spikes = zeros(n_frames,1);
        display('getting frames')
        tic
        for f = 1:n_frames
            
            if SU
                [Spike_Timing index numtrials] = getTrials(frameEpocs{blocknum},times, f, frame_duration);
                eps = frameEpocs{blocknum};
            else
                [Spike_Timing index numtrials] = getTrialsSU(data.frameEpocs,data.MUspikeT{cell_n}, f, frame_duration);
                eps = data.frameEpocs;
            end
            n_spikes(f) = length(Spike_Timing);
            ntrials(f) = numtrials;
        end
        toc
        
        n_reps = min(ntrials);
        
        
        %%% evaluate these before normalizing to rate
        N=nansum(n_spikes)
        duration_wn(cell_n)=max(times);
        sta_N(cell_n)=N;
        fano(cell_n) = var(n_spikes)/mean(n_spikes);
        
        if N==0
            break
        end
        n_spikes = n_spikes./(frame_duration*ntrials');
        n_spikes(isnan(n_spikes))=0;
        
        
        if movietype == cm_noise
            wnfig=figure;
            subplot(2,2,1)
            cspikes=condenseData(n_spikes,framerate);
            bar(cspikes);
            xlim([0 length(cspikes)]) ;
            xlabel('sec');ylabel('sp/sec')
        end
        
        cyc_frames = round(contrast_period/frame_duration);
        
        if contrast_modulated
            %     for f= 1:n_frames;
            %         frm = moviedata(:,:,f);
            %         c(f) = double(std(frm(:)));
            %     end
            f = 1:n_frames;
            c = double(0.5- 0.5*cos(2*pi*f/cyc_frames));  %%% define contrast
            
            f_null=find(c==0);
            spont = sum(n_spikes(f_null).*ntrials(f_null)')/sum(ntrials(f_null));
            wn(cell_n,stim_eye).spont=spont;
            
            nf = zeros(cyc_frames,1);
            contrastdata = zeros(cyc_frames,1);
            for i = 1:n_frames;
                contrastdata(mod(i,cyc_frames)+1)=contrastdata(mod(i,cyc_frames)+1)+double(c(i));
                nf(mod(i,cyc_frames)+1) = nf(mod(i,cyc_frames)+1)+1;
            end
            contrastdata = contrastdata./nf;
            
            fftsignal = (fft(n_spikes));
            fftsignal = fftsignal(1:round(n_frames/2));
            subplot(2,2,2)
            loglog(abs(fftsignal));
            
            df=(1/n_frames);
            contrastfreq = 1/cyc_frames;
            fft_chan = round(n_frames/cyc_frames)+1   %%%contrastfreq/df
            framerate = round(1/frame_duration);
            
            
            
            wn(cell_n,stim_eye).responsiveness = 2*abs(fftsignal(fft_chan))/(abs(fftsignal(1)));  %%double to normalize FFT relative to DC
            wn(cell_n,stim_eye).phase = angle(fftsignal(fft_chan));
            phaseangle = angle(fftsignal(fft_chan))
            cycledata = zeros(cyc_frames,1);
            for i = 1:n_frames;
                cycledata(mod(i,cyc_frames)+1)=cycledata(mod(i,cyc_frames)+1)+n_spikes(i);
            end
            cycledata = cycledata./(nf);
            
            %cycledata = conv(cycledata,ones(1,30))/10;
            cycledata = condenseData(cycledata,framerate/2);
            contrastdata = condenseData(contrastdata,framerate/2);
            
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
            %figure
            subplot(2,2,3)
            hold on
            %%% plot average of up and down on contrast response
            meancontrast = 0.5*(contrastdata(1:10)+contrastdata(20:-1:11));
            meandata = 0.5*(cycledata(1:10)+cycledata(20:-1:11));
            plot(meancontrast/max(meancontrast),meandata,'g');
            title(sprintf('unit %d %d',channel_no, clust_no))
            xlabel('contrast');
            ylabel('sp / sec');
            plot(contrastdata/max(contrastdata),cycledata);
            
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
        %     saveas(gcf,fullfile (pname,sprintf('hist%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
        
        
        %%% number of time points to calculate STA at
        
        
        clear sp
        
        display('calculating sta')
        tic
        for lag = 1:sta_length
            sp(:,lag) = n_spikes(lag:end-sta_length+lag);
        end
        if movietype==cm_noise
            st = mov_squeeze*sp;
            sta_t = reshape(st,size(moviedata,1),size(moviedata,2),sta_length)/sum(n_spikes);
        end
        N= sum(n_spikes);
        toc
        
        if movietype == mv_noise | movietype==fl_noise
            sp=single(sp);
            sta_neg =  reshape(mov_neg*sp,size(moviedata,1),size(moviedata,2),sta_length)/N;
            sta_zero=reshape(mov_zero*sp,size(moviedata,1),size(moviedata,2),sta_length)/N;
            sta_pos = reshape(mov_pos*sp,size(moviedata,1),size(moviedata,2),sta_length)/N;
            sta_t = sta_pos-sta_neg;
            sta_diff = sqrt((sta_neg-mov_neg_avg).^2 + (sta_zero-mov_zero_avg).^2 +(sta_pos-mov_pos_avg).^2);
        end
        
        if pos_neg
            moviepos = (moviedata-127)/127;
            moviepos(moviepos<0)=0;
            moviepos = reshape(moviepos(:,:,1:end-sta_length),size(moviedata,1)*size(moviedata,2),size(moviedata,3)-sta_length);
            st = moviepos*sp;
            sta_pos = reshape(st,size(moviedata,1),size(moviedata,2),sta_length)/N;
            
            movieneg = (moviedata-127)/127;
            movieneg(movieneg>0)=0;
            movieneg = reshape(movieneg(:,:,1:end-sta_length),size(moviedata,1)*size(moviedata,2),size(moviedata,3)-sta_length);
            st = movieneg*sp;
            sta_neg = reshape(st,size(moviedata,1),size(moviedata,2),sta_length)/N;
            
        end
        
        if movietype ==cm_noise
            color_range = [-64 64];
            movavg=127;
        else
            color_range = [-0.25 0.25];
            movavg=0;
        end
        %%% plot STA at each time point
        stafig=figure
        for t = 1:16
            subplot(4, 4, t);
            
            imagesc(sta_t(:,:,t)'-movavg' ,color_range);
            axis equal
            axis tight
            set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])
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
        drawnow;
        
        if movietype == mv_noise | movietype==fl_noise
            figure
            for t = 1:16
                subplot(4, 4, t);
                imagesc(sta_diff(:,:,t)' ,[0 .15]);
                axis equal
                axis tight
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
                
                if t==1
                    title_text = sprintf('%s',Tank_Name);
                    text(0,-10,title_text,'FontSize',8);
                end
                if t==2
                    title_text = sprintf('ch%d c%d',channel_no,clust_no);
                    text(0,-10,title_text,'FontSize',8);
                end
            end
            drawnow;
        end
        
        if movietype~=cm_noise
            
            clear d
            
            %[m ind] = max(abs(sta_t(:)-127));
            
            [m ind] = max(sta_diff(:));
            m
            [x y t] = ind2sub(size(sta_t),ind)
            t_lag = t-1;
            
            if movietype == mv_noise
                d=zeros(4,size(moviedata,3)-t_lag-1);
            else
                d=zeros(2,size(moviedata,3)-t_lag-1);
            end
            size(d)
            
            clear p
            p{1} = [-1 1];
            if movietype==fl_noise
                p{2} = [1 2 4 8 16 255];
            else
                p{2} = [2 4 8];
            end
            %%% contrast and size at point x,y over course of movie
            d(1,:) = squeeze(round((double(moviedata(x,y,1:end-t_lag-1))-127)/127));
            d(2,:) = squeeze(sz_mov(x,y,1:end-t_lag-1));
            hist_all = zeros(length(p{1}),length(p{2}));
            if movietype==mv_noise
                p{3} = [10 20 40 80 160];
                p{4} =1:8;
                d(3,:) = squeeze(sp_mov(x,y,1:end-t_lag-1));
                d(4,:) = th_mov(x,y,1:end-t_lag-1);
                hist_all = zeros(length(p{1}),length(p{2}),length(p{3}),length(p{4}));
            end
            
            size(d)
            display('computing histogram');
            tic
            
            n_all = hist_all;
            ndim = length(p);
            if movietype== fl_noise
                for c = 1:length(p{1})
                    c
                    for sz = 1:length(p{2})
                        f = find(d(1,:)==p{1}(c) & d(2,:)==p{2}(sz));
                        n_all(c,sz)=length(f);
                        if length(f)>0
                            hist_all(c,sz)=sum(n_spikes(f+t_lag).*ntrials(f+t_lag)')/sum(ntrials(f+t_lag));
                        else
                            hist_all(c,sz)=nan;
                        end
                    end
                end
            end
            
            if movietype == mv_noise
                for c = 1:length(p{1})
                    c
                    for sz = 1:length(p{2})
                        for sp = 1:length(p{3})
                            for th = 1:length(p{4})
                                hit = (d(1,:)==p{1}(c) & d(2,:)==p{2}(sz) & d(3,:)==p{3}(sp) & d(4,:)==p{4}(th));
                                n_unique(c,sz,sp,th) = sum(diff(hit)>0);
                                f= find(hit);
                                n_all(c,sz,sp,th)=length(f);
                                if length(f)>0
                                    hist_all(c,sz,sp,th)=sum(n_spikes(f+t_lag).*ntrials(f+t_lag)')/sum(ntrials(f+t_lag));
                                else
                                    hist_all(c,sz,sp,th)=nan;
                                end
                            end
                        end
                    end
                end
            end
            
            squeeze(hist_all(1,:,:,:))
            squeeze(hist_all(2,:,:,:))
            toc
            
            f_null = find(d(1,:)==0);
            spont = sum(n_spikes(f_null+t_lag).*ntrials(f_null+t_lag)')/sum(ntrials(f_null+t_lag));
            
            if movietype==mv_noise
                figure
                h = (squeeze(nanmean(hist_all,4)));
                h(isnan(h))=-0.1*max(h(:));
                for i= 1:2
                    subplot(1,2,i)
                    
                    imagesc(squeeze(h(i,:,:)),[min(h(:)) max(h(:))]);
                    axis equal
                    axis square
                end
            end
            
            if movietype == fl_noise
                sz_tune= hist_all;
            elseif movietype == mv_noise
                sz_tune = nanmean(nanmean(hist_all,4),3);
            end
            
            tuningfig=figure
            subplot(2,2,1);
            plot(sz_tune');
            hold on
            plot([1 length(sz_tune)],[spont spont],'r')
            title(sprintf('ch%d c%d',channel_no,clust_no))
            yl = ylim;
            axis([1 length(sz_tune) 0 yl(2)])
            xlabel('size')
            set(gca,'Xtick',1:length(p{2}));
            set(gca,'Xticklabel',p{2});
            
            if movietype == mv_noise
                
                
                sp_tune = squeeze(nanmean(nanmean(hist_all,4),2));
                subplot(2,2,2);
                plot(sp_tune');
                hold on
                plot([1 5],[spont spont],'r')
                yl=ylim;
                axis([1 5 0 yl(2)])
                xlabel('speed');
                title(printf('ch%d c%d',channel_no,clust_no))
                set(gca,'Xtick',1:length(p{3}));
                set(gca,'Xticklabel',p{3});
                
                th_tune = squeeze(nanmean(nanmean(hist_all,3),2));
                
                subplot(2,2,3);
                plot(th_tune');
                hold on
                plot([1 8],[spont spont],'r')
                yl = ylim;
                axis([1 8 0 yl(2)])
                xlabel('theta');
                title(printf('ch%d c%d',channel_no,clust_no))
                
                
            end
            subplot(2,2,4);
            imagesc(sta_t(:,:,t_lag+1)',[-0.25 0.25]);hold on
            plot(x,y,'o');
            xlabel(sprintf('t = %d',t_lag+1));
            
            
            %%%%% histrograms relative to onset/offset
            timefig = figure;
            sizes = [1 2 4 8 16 255];
            for rep = 1:2
                for sz = 1:length(sizes);
                    n=0;
                    ontime=0;
                    tseries = single(squeeze(moviedata(x,y,:)));
                    size_series = single(squeeze(sz_mov(x,y,:)));
                    if rep ==1
                        if movietype==mv_noise
                            onset = (tseries(1:end-1)==127 & diff(tseries)<0 &size_series(2:end)==sizes(sz));
                        else
                            onset = (tseries(1:end-2)==127 & tseries(2:end-1)==0 & tseries(3:end)==127 &size_series(2:end-1)==sizes(sz));
                        end
                        c = 'b';
                    elseif rep ==2
                        if movietype==mv_noise
                            onset = (tseries(1:end-1)==127 & diff(tseries)>0 &size_series(2:end)==sizes(sz));
                        else
                            onset = (tseries(1:end-2)==127 & tseries(2:end-1)==255 & tseries(3:end)==127 &size_series(2:end-1)==sizes(sz));
                        end
                        c= 'r';
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
                        dt=0.05;
                        histbins = 0:dt:0.8;
                    elseif movietype==mv_noise
                        dt = 0.1;
                        histbins = -1:dt:2;
                    end
                    
                    h=0;
                    subplot(2,1,1)
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
                        clear onset_hist
                        
                        subplot(2,1,2)
                        hold on
                        plot(histbins(1:end-1)+dt/2,h(1:end-1)/(n*dt),c)
                        onset_hist(rep,sz,:) = h(1:end-1)/(n*dt);
                        onset_bins=histbins;
                    else
                        onset_hist=NaN
                        onset_bins=histbins
                    end
                end
            end
            subplot(2,1,2)
            legend('gray->on','gray->off')
            title(sprintf('ch%d c%d',channel_no,clust_no))
            
            if movietype ==fl_noise
                plot([0.25 0.5], [0.5 0.5],'g','LineWidth',8)
            end
            
            
            
        end
        
        
        
        nx = size(sta_t,1);
        ny = size(sta_t,2);
        nt = size(sta_t,3);
        
        
        
        if compute_svd
            sta_col =reshape(sta_t,nx*ny,nt);
            [u s v] = svd(sta_col-mean(sta_t(:)));
            
            svdfig=figure
            
            for i = 1:3
                subplot(2,3,i)
                range = max(abs(min(u(:,i))),abs(max(u(:,i))));
                imagesc(reshape(u(:,i),nx,ny)'*sign(v(4,i)),1.2*[-range range]);
                svdt(i,:,:) = reshape(u(:,i),nx,ny)'*sign(v(4,i));
                v(:,i) = v(:,i)*sign(v(4,i));
                axis equal; axis tight;
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
                subplot(2,3,i+3)
                plot(v(:,i));
                xlim([0 size(v,1)])
            end
            
            %figure
            %plot(diag(s));
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
            
            [x y t_lag] = ind2sub(size(sta_t),ind)
            figure(wnfig)
            subplot(2,2,4)
            imagesc(fftshift(abs(fft2(sta_t(:,:,t_lag)'-128))));
            title(titlestr);
            axis equal;
            set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])
            axis tight
            
            
            
            %saveas(gcf,fullfile(noisepname,sprintf('fft_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
            
            %             h = zeros(size(n_spikes));
            %             m=0;
            %             tic
            %             sta=(sta_t(xrange,yrange,t_lag)-movavg(xrange,yrange))/128;
            %
            %
            %             sta = sta/sqrt(sum(sum(sta.^2)));
            %
            %             %%% calculate transfer function
            %             for t = (t_lag+1):n_frames
            %                 if n_spikes(t)>=0
            %                     m = (moviedata(xrange,yrange,t-t_lag+1)-movavg(xrange,yrange))/(128*sqrt(dx*dy));
            %                     h(t) = sta(:)'*m(:);
            %                 end
            %             end
            %
            %             %      figure
            %             %      plot(h(t_lag+1:n_frames),n_spikes(t_lag+1:n_frames),'.');
            %
            %             hist_int =0.05;
            %             h_round = round(h(t_lag+1:n_frames)/hist_int);
            %             n_sp = n_spikes(t_lag+1:n_frames);
            %             h_min = min(h_round);
            %             h_range = h_min:max(h_round)-1;
            %             n_mean=0; n_std=0; n_samp=0;
            %             for i = h_range;
            %                 use = find(h_round ==i);
            %                 n_mean(i-h_min+1) = mean(n_sp(use));
            %                 n_std(i-h_min+1) = std(n_sp(use));
            %                 n_samp(i-h_min+1) = size(use,1);
            %             end
            %             %             if movietype == cm_noise
            %             %                 figure
            %             %                 errorbar(h_range*hist_int,n_mean,n_std./sqrt(n_samp));
            %             %                 hold on;
            %             %                 plot(h_range*hist_int,n_samp/500,'g');
            %             %             end
            %
            %             transfer_function(cell_n,1:size(n_mean,2)) = n_mean;
        end
        
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
        
        if movietype == cm_noise
            
            saveas(wnfig,fullfile(noisepname,sprintf('wndata_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
            saveas(stafig,fullfile(noisepname,sprintf('wn_sta_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
            saveas(svdfig,fullfile(noisepname,sprintf('wn_svd_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
            
            figure(wnfig);
            title(sprintf('ch=%d cl=%d',channel_no,clust_no));
            set(gcf, 'PaperPositionMode', 'auto');
            print('-dpsc',psfilename,'-append');
            
            figure(stafig);
            set(gcf, 'PaperPositionMode', 'auto');
            print('-dpsc',psfilename,'-append');
            
            figure(svdfig);
            set(gcf, 'PaperPositionMode', 'auto');
            print('-dpsc',psfilename,'-append');
            
            wn(cell_n,stim_eye).crf=cycledata;
            wn(cell_n,stim_eye).N=N;
            wn(cell_n,stim_eye).svd_xy = svdt;
            wn(cell_n,stim_eye).svd_t = v;
            wn(cell_n,stim_eye).sta=sta_t;
            %wn(cell_n,stim_eye).spont=spont;
            
            close all
        elseif movietype == fl_noise
            
            saveas(timefig,fullfile(noisepname,sprintf('flash_time_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
            saveas(tuningfig,fullfile(noisepname,sprintf('flash_tuning_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
            %saveas(tuningfig,fullfile(noisepname,sprintf('flash_tuning_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
            
            figure(timefig)
            title(sprintf('ch=%d cl=%d',channel_no,clust_no));
            set(gcf, 'PaperPositionMode', 'auto');
            print('-dpsc',psfilename,'-append');
            
            figure(tuningfig)
            
            set(gcf, 'PaperPositionMode', 'auto');
            print('-dpsc',psfilename,'-append');
            
            fl(cell_n).N=N;
            fl(cell_n).hist_all=hist_all;
            fl(cell_n).n_all = n_all;
            fl(cell_n).sta_diff=sta_diff;
            fl(cell_n).sta=sta_t;
            fl(cell_n).lag=t_lag;
            fl(cell_n).spont=spont;
            fl(cell_n).onset_hist=onset_hist;
            fl(cell_n).onset_bins=onset_bins;
            fl(cell_n).sta_pos=[x y];
            
   if ~isnan(onset_hist);
       
               OnOffFlash;
     fl(cell_n).sustainBias=sustain;
     fl(cell_n).onoffbias = onoffbias;
     fl(cell_n).RFzcore = rfz;
     fl(cell_n).rf = rf;
     fl(cell_n).resps = resps; %%% timecourse averaged over sizes
     fl(cell_n).flash_resp = flash_resp;  %%% mean response for on/off and size
   else
     fl(cell_n).sustainBias=NaN;
     fl(cell_n).onoffbias = NaN;
     fl(cell_n).RFzcore = NaN;
     fl(cell_n).rf = NaN;
     fl(cell_n).resps = NaN; %%% timecourse averaged over sizes
     fl(cell_n).flash_resp = NaN;
     
   end
        elseif movietype==mv_noise
            
            saveas(timefig,fullfile(noisepname,sprintf('move_time_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
            saveas(tuningfig,fullfile(noisepname,sprintf('move_tuning_%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
            
            figure(timefig)
            title(sprintf('ch=%d cl=%d',channel_no,clust_no));
            set(gcf, 'PaperPositionMode', 'auto');
            print('-dpsc',psfilename,'-append');
            
            
            figure(tuningfig)
            
            set(gcf, 'PaperPositionMode', 'auto');
            print('-dpsc',psfilename,'-append');
            
            mv(cell_n).N=N;
            mv(cell_n).hist_all=hist_all;
            mv(cell_n).n_unique = n_unique;
            mv(cell_n).sta_diff=sta_diff;
            mv(cell_n).sta=sta_t;
            mv(cell_n).lag=t_lag;
            mv(cell_n).spont=spont;
            mv(cell_n).onset_hist=onset_hist;
            mv(cell_n).onset_bins=onset_bins;
            mv(cell_n).sta_pos=[x y];
        end
        
        
    end

close all
%         OnOffFlash
end  %%%cell


delete(gcp('nocreate'))

ps2pdf('psfile', psfilename, 'pdffile', [psfilename(1:(end-2)) 'pdf']);
delete(psfilename);


% parpool close
% 

if movietype==cm_noise
    if ~exist ('sess','var')
    sess = input('post doi? 0/1 : ')+1;
    end
    
    wn(cell_n,stim_eye).degperpix=degperpix;
    load(afile,'wnblocks');
    wnblocks{sess}=wn;
    save(afile,'wnblocks','-append');
    if sess==2
        wn_post=wn;
        save(afile,'wn_post','-append')
    elseif sess==1;
        save(afile,'wn','-append')
    end
elseif movietype==fl_noise
    fl(cell_n).degperpix=degperpix;
    save(afile,'fl','-append')
elseif movietype==mv_noise
    mv(cell_n).degperpix=degperpix;
    save(afile,'mv','-append')
end





