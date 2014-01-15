%function drift_orientfreq_laser
% Matlab codes for reading from TTank for drifting or counterphase gratings
% plots histgrams and rasters and performs fourier analysis
%function gratfreq_cluster
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03

%%% read in cluster data, then connect to the tank and read the block



% Tank_Name='06222011_mlr_stim_awake'
% Block_Name='wn3_72hz'

nchan = input('# chans : ');
chans = 1:4:nchan;

SU = input('multiunit (0) or single-unit (1) : ');
if SU
    [fname, pname] = uigetfile('*.mat','cluster data');
    load(fullfile(pname,fname));
    for i =1:length(Block_Name);
        sprintf('%d : %s ',i,Block_Name{i})
    end
    block = input('which block to analyze ? ');
    Block_Name = Block_Name{block}
    [afname, apname] = uigetfile('*.mat','analysis data');
    noisepname = apname;
    afile = fullfile(apname,afname);
    load(afile);
    use_afile=1;
    cells
else
    pname = uigetdir('C:\data\TDT tanks','block data')
delims = strfind(pname,'\');
selected_path = pname(1 :delims(length(delims))-1)
Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
Block_Name = pname(delims(length(delims))+1 :length(pname))
end

laser = input('movement (0) or laser (1) : ');

if laser
    flags = struct('laserOn',1,'MUspike',1,'visStim',1)
    tdtData= getTDTdata(Tank_Name, Block_Name, chans, flags);
    tsamp  = tdtData.laserT;
    vsmooth = tdtData.laserTTL;
else
    flags = struct('mouseOn',1,'MUspike',1,'visStim',1)
    tdtData= getTDTdata(Tank_Name, Block_Name, chans, flags);
    tsamp = tdtData.mouseT;
    vsmooth = tdtData.mouseV;
end




thresh_velocity = 1;
figure
plot(tsamp,vsmooth);


plot_duration=2.5; %in second

hist_int = 0.05;
hist_range=[0:hist_int:plot_duration];
axis_range=[0 plot_duration 0 25];
max_events=50000;

stim_duration =1.5;
wait_duration = .5;  %duration after each stimulus (for calculating spontaneous rate)
blank_interval = 0.1; %% length of time after stimulus not to use in calculating spontaneous rate
fft_int = .05;  %%% size of bins to use for fourier analyis (in sec)
tempfreq = 2;    %%% temporal frequency of stimulus
blank_stim = 1; %%% is there an extra stimulus to measure spontaneous
full_field = 1;  %%% is there a full-field flicker?
tetrode_linear=0;

if tetrode_linear
    ch_map = [14 8 10 4 13 7 9 3 11 1 15 5 12 2 16 6];
else
    ch_map = 1:16;
end


% set number of conditions and display setup (generally rows = orientation, columns = frequency)
n_rows=12;
n_col=6;
orients = [0 45 90 135 180 225 270 315];
orients = 0:30:330;
spatfreqs = [.01 .02 .04 .08 .16 .32];
%



n_cond=n_rows*n_col;
if blank_stim
    n_cond=n_cond+1;
end
if full_field
    n_cond=n_cond+1;
end

%printfig = input('print ? ');
printfig=0;

if SU
    cell_range = 1:size(cells,1)
else
    cell_range=1:4:nchan;
end
for cell_n = cell_range;
    % for cell_n=9:9
    cell_n
    if SU
        channel_no = cells(cell_n,1)
        clust_no = cells(cell_n,2)
        channel_times =spikeT{cell_n} - (block-1)*10^5;
        times = channel_times(channel_times>0 & channel_times<10^5);
        hist_fig = figure('Name',sprintf('unit %d %d',channel_no,clust_no))
    else
        hist_fig = figure('Name',sprintf('channel %d',cell_n))
        channel_no = cell_n;
        clust_no = [];
    end
    
    for rep =1:2
        hist_fig = figure;
        rast_fig = figure;
        % fft_fig = figure;
        %         spont_full= figure;
        spont_full_rast=figure;
        for cond =0:n_cond-1
                  if SU
                [Spike_Timing index numtrials Epocs_TS] = getTrialsSU(stimEpocs{block},times, cond+1, stim_duration);
            
            else
               [Spike_Timing index numtrials Epocs_TS] = getTrialsSU(tdtData.stimEpocs,tdtData.MUspikeT{cell_n}, cond+1, stim_duration);
                  end
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% movement specific
            clear trial_velocity newtrial
            for i = 1:numtrials;
                trial_velocity(i) = mean(vsmooth(find(tsamp>Epocs_TS(i) & tsamp<Epocs_TS(i)+stim_duration)));
            end
            trial_velocity
            if rep==1
                usedtrial = trial_velocity<thresh_velocity;
            else
                usedtrial = trial_velocity>thresh_velocity;
            end
        
            trials = find(usedtrial)                   
            if isempty(trials)
                [m trials] = min(trial_velocity);
                trials
                usedtrial(trials)=1;
            end
        
            for i = 1:numtrials;
                if isempty(find(trials==i))
                    newtrial(i)=0;
                else
                    newtrial(i) = find(trials==i);
                end
            end
       
            numtrials = sum(usedtrial);
         %%% only keep spikes that are in the desired cluster

            Spike_Timing = Spike_Timing(find(usedtrial(index)));
            index=index(usedtrial(index));
            index = newtrial(index);  %%% get rid of unused trials
            
            %%%% end movememnt specific
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            title_text=['Orientation: ' num2str(cond*45)];
            
            %raster plot
            if cond<n_rows*n_col;
                figure(rast_fig);
                subplot(n_rows,n_col,cond+1); hold on; set(gca, 'yDir','reverse');
                
                %  plot (Spike_Timing, index, '.k', 'MarkerSize',4);
                plot ([Spike_Timing; Spike_Timing], [index-0.25;index+0.25], 'k', 'MarkerSize',4);
                
                axis([0 plot_duration 0 numtrials+1]);
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
                if cond==0
                    title_text = sprintf('%s',Tank_Name);
                    text(0,3,title_text,'FontSize',8);
                end
                if cond==1
                    title_text = sprintf('ch%d c%d',channel_no,clust_no);
                    text(0,3,title_text,'FontSize',8);
                end
            end
            %% histograms
            rate_hist = hist(Spike_Timing, hist_range)/(hist_int*numtrials);
            if cond<n_rows*n_col
                figure(hist_fig);
                subplot(n_rows,n_col,cond+1);
                bar(hist_range, rate_hist);  hold on;
                axis(axis_range);
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
                if cond==0
                    title_text = sprintf('%s',Tank_Name);
                    text(0,20,title_text,'FontSize',8);
                end
                if cond==1
                    title_text = sprintf('ch%d c%d',channel_no,clust_no);
                    text(0,20,title_text,'FontSize',8);
                end
            end
            
            if cond>=n_rows*n_col & full_field%%% blank frame or full-field, but only plot if full-field is on
                %                 figure(spont_full)
                %                 subplot(2,1,cond-n_rows*n_col +1);
                %                 bar(hist_range, rate_hist);  hold on;
                %                 axis(axis_range);
                %                 set(gca,'XTickLabel',[])
                %                 set(gca,'YTickLabel',[])
                %                 title(title_text);
                %
                figure(spont_full_rast);
                
                subplot(2,1,cond-n_rows*n_col +1);
                hold on; set(gca, 'yDir','reverse');
                plot (Spike_Timing, index, '.k', 'MarkerSize',8);
                axis([0 plot_duration 0 numtrials+1]);
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
                title(title_text);
            end
            
            
            %%fourier analysis
            fft_range = 0:fft_int:stim_duration-fft_int;
            fft_range = fft_range + (fft_int/2);  %%% to center bins
            
            nyq = 0.5/fft_int; %%% nyquist interval
            freq_int = nyq / (0.5*stim_duration/fft_int);
            
            
            fft_data=abs(fft(hist(Spike_Timing(Spike_Timing<stim_duration),fft_range)/(fft_int*numtrials)));
            fft_data = fft_data/size(fft_data,2);
            %fft_data(1) = fft_data(1)/2;
            %             if cond<n_rows*n_col;
            %                 figure(fft_fig);
            %                 subplot(n_rows,n_col,cond+1);
            %                 plot(fft_data(1:round((size(fft_data,2)+1)/2)));
            %                 axis([0 16 0 8])
            %             end
            
            cond
            numtrials
            R(cell_n,cond+1) = sum(Spike_Timing<=stim_duration)/(stim_duration*numtrials);
            spont(cell_n,cond+1) = sum((Spike_Timing>(stim_duration+blank_interval))&(Spike_Timing<(stim_duration+wait_duration)))/(numtrials*(wait_duration-blank_interval));
            
            F0(cell_n,cond+1) = fft_data(1);
            F1(cell_n,cond+1) = 2*fft_data(1+tempfreq/freq_int);  %%% double to count both pos & neg frequency
            F2(cell_n,cond+1) = 2*fft_data(1+2*tempfreq/freq_int);
            
            clear spikeR
            for i = 1:numtrials;
                spikeR(i) = sum(Spike_Timing<=stim_duration & index==i)/(stim_duration);
            end
            R_var(cell_n,cond+1) = var(spikeR);
            R_err(cell_n,cond+1) = std(spikeR)/sqrt(numtrials -1);
            cond
            mean(spikeR)
            var(spikeR)
            std(spikeR)/sqrt(numtrials -1)
            
            spikeR
        end  %% cond
        
        title_text = sprintf('channel %d cluster %d',channel_no, clust_no);
        
        if blank_stim
            spont_avg = R(cell_n,n_rows*n_col+1);
            both_drift_spont(cell_n,rep)=R(cell_n,n_rows*n_col+1);
            % spont_avg = 0.5*( R(cell_n,n_rows*n_col+1) +  mean(spont(cell_n,:)))
            
        else
            spont_avg = mean(spont(cell_n,:))  %%% if no blank frame, average over inter-stimulus interval of all conditions
        end
        
        
        
        
        plotcolor = 'bgrcmykbgr';
        
        tuning_fig = figure;
        
        
        %for f= 1:n_col;
        for f = 1:n_col
            plot(R(cell_n,f:n_col:f+n_col*(n_rows-1))-spont_avg,plotcolor(f));
            hold on;
        end
        title(title_text);
        legend('.01 cpd','.02cpd','.04cpd','.08cpd','.16cpd','.32cpd')
        %set(gca,'XTickLabel',['0' '45' '90' '135' '180' '225' '270' '315']);
        
        
        %% calculate tuning parameters
        
        orientfreq = reshape(R(cell_n,1:n_col*n_rows),n_col,n_rows)'-spont_avg;
        orient_tuning_all = mean(orientfreq,2);
        %  figure
        %  plot(orient_tuning_all);
        if max(orient_tuning_all)>abs(min(orient_tuning_all))
            [max_resp pref_orient(cell_n)] = max(orient_tuning_all);
        else
            [max_resp pref_orient(cell_n)] = min(orient_tuning_all);
        end
        freq_tuning(cell_n,:) = orientfreq(pref_orient(cell_n),:);
        figure
        subplot(2,1,1);
        plot(freq_tuning(cell_n,:));
        if max(freq_tuning)>abs(min(freq_tuning))
            [max_resp pref_freq(cell_n)] = max(freq_tuning(cell_n,:));
        else
            [max_resp pref_freq(cell_n)] = min(freq_tuning(cell_n,:));
        end
        orient_tuning(cell_n,:) = orientfreq(:,pref_freq(cell_n))';
        subplot(2,1,2);
        plot(orient_tuning(cell_n,:));
        % f1f0(cell_n) = F1(cell_n,(pref_orient-1)*n_col +pref_freq)/ F0(cell_n,(pref_orient-1)*n_col +pref_freq)
        
        
        
        %%% reshape, and use average across SF to determine preferred angle
        n_orient = n_rows;
        orientfreq(:,2:n_col+1) = reshape(R(cell_n,1:n_col*n_rows),n_col,n_rows)'-spont_avg;
        orientfreq(:,1) = R(cell_n,n_col*n_rows+2)-spont_avg;  %%% full-field (w=0)
        
        n_w = n_col+1;
        orient_tuning_mn = mean(orientfreq,2)';
        
        theta_ind = orients*pi/180;
        [theta_pref OSI A1 A2 w B yfit] = fit_tuningcurve(orient_tuning_mn,theta_ind);
        
        %       figure
        %      plot(orient_tuning_mn);
        %        hold on
        %        plot(yfit,'g')
        %          title(title_text);
        
        sf = [0:n_col];
        %%% find preferred SF at optimal theta
        [w theta] = meshgrid(sf,theta_ind);
        
        sf_spacing=0.25;
        interp_sfs =0:sf_spacing:6;
        theta_pref_col = ones(size(interp_sfs));
        
        theta_pref_col(:) = theta_pref;
        %%% add theta =2*pi to orientfreq and theta, so we don't interpolate past
        %%% the end for theta_pref>330 deg
        orientfreq_w = orientfreq;
        orientfreq_w(n_rows+1,:) = orientfreq_w(1,:);
        theta_ind_all = theta_ind;
        theta_ind_all(length(theta_ind)+1)=2*pi;
        [w_all theta_all] = meshgrid(sf,theta_ind_all);
        w_tuning_curve = interp2(w_all,theta_all,orientfreq_w,interp_sfs,theta_pref_col,'bicubic');
        
        theta_pref_row = round(theta_pref/(pi/6))+1;
        if theta_pref_row > length(orients)
            theta_pref_row =1;
        end
        theta_pref_row
        w_tuning_round = interp1(sf,orientfreq(theta_pref_row,:),interp_sfs)   ;
        
        %     figure
        %     plot(w_tuning_curve);
        %     hold on
        %     plot(w_tuning_round,'r');
        %%% if can't interpolate tuning curve, use rounded theta
        if isnan(w_tuning_curve)
            w_tuning_curve= w_tuning_round;
        end
        
        drift_sftuning(cell_n,:,rep)=w_tuning_curve;
        %         [amp w_pref] = max(w_tuning_curve);
        %         if w_pref>(1/sf_spacing);
        %             [amp w_pref w_bw B yfit] = fit_1d_gaussian(w_tuning_curve);
        % %             hold on
        % %             plot(yfit,'g');
        % %             plot(1:4:25,w_tuning_curve(1:4:25),'o');
        %
        %         else
        %             w_bw=0;
        %         end
        
        [wpref_dog(cell_n,rep) wbw_dog(cell_n,rep) peak_value(cell_n,rep)] = dogfit(squeeze(drift_sftuning(cell_n,:,rep)));
        wpref_dog(cell_n,rep)
        w_pref = log2(wpref_dog(cell_n,rep)/.01)+1
        if w_pref<=0
            w_pref=0.5;
        end
        if w_pref>6
            w_pref=6;
        end
        w_pref   %%% this is the index (i.e. column) of preferred w
        w_bw = wbw_dog(cell_n,rep)
        
        title(title_text);
        
        
        %%% fit orientation tuning curve at optimal SF
        w_pref_col = ones(size(1:n_orient)).*w_pref;
        theta_tuning = interp2(w,theta,orientfreq,w_pref_col,theta_ind);
        
        [driftorient_theta(cell_n,rep) driftorient_osi(cell_n,rep) driftorient_A1(cell_n,rep) driftorient_A2(cell_n,rep) ...
            driftorient_thetawidth(cell_n,rep) driftorient_B(cell_n,rep) driftorient_null(cell_n,rep) yfit]= fit_tuningcurve(theta_tuning,theta_ind);
        
        i=sqrt(-1);
        mu = (sum(theta_tuning.*exp(i*theta_ind)))/sum(theta_tuning);
        if isnan(mu);
            mu=0;
        end
        
        driftorient_dsi(cell_n) = abs(mu);
        
        figure
        plot(theta_ind,theta_tuning);
        hold on
        plot(theta_ind,yfit(1:length(theta_ind)),'g');
        title(title_text);
        
        if w_pref <1
            driftorient_wpref(cell_n)=0;
        else
            driftorient_wpref(cell_n) = .01*(2^(w_pref-1))
        end
        driftorient_bw(cell_n) = w_bw;
        
        drift_thetatuning(cell_n,:,rep) = theta_tuning;
        drift_sftuning(cell_n,:,rep) = w_tuning_curve;
        drift_sftuning_round(cell_n,:) = w_tuning_round;
        
        
        
        %%% interpolate F0, F1 values
        theta_pref = driftorient_theta(cell_n);
        %% to avoid interpolating past end of theta
        %%%(for angles between 2pi-pi/8 and  and 2pi
        if theta_pref > max(max(theta))
            theta_pref = max(max(theta));
        end
        
        F0all(:,2:n_col+1) = reshape(F0(cell_n,1:n_col*n_rows),n_col,n_rows)';
        F0all(:,1) = F0(cell_n,n_col*n_rows+2);  %%% full-field (w=0);
        driftorient_F0(cell_n,rep) = interp2(w,theta,F0all,round(w_pref), theta_pref );
        
        F1all(:,2:n_col+1) = reshape(F1(cell_n,1:n_col*n_rows),n_col,n_rows)';
        F1all(:,1) = F1(cell_n,n_col*n_rows+2);  %%% full-field (w=0)
        driftorient_F1(cell_n,rep) = interp2(w,theta,F1all,round(w_pref), theta_pref );
        
        drift_spont(cell_n,rep) = R(cell_n,n_cond-1);
        drift_orientfreq_all(cell_n,:,:)=orientfreq;
        
        [drift_peakR(cell_n) maxcond] = max(R(cell_n,1:n_cond-2));
        drift_peakvar(cell_n) = R_var(cell_n,maxcond);
        drift_peakerr(cell_n) = R_err(cell_n,maxcond);
        
        title_text = sprintf('%s ch%d cl%d',Tank_Name, channel_no,clust_no);
        
        if printfig
            print(hist_fig);
        end
        %         figure(fft_fig);
        %         subplot(12,8,1);
        %         text(0,10,title_text);
        %         figure(spont_full);
        %         title(title_text);
        %         figure(spont_full_rast);
        %         title(title_text);
        % close(spont_full)
        % close(spont_full_rast)
        %
        xlabel('secs');
        % saveas(tuning_fig,fullfile(pname,sprintf('grattuning_move%d%s_%d_%d',rep,Block_Name,channel_no,clust_no)),'fig')
%         saveas(rast_fig,fullfile(apname,sprintf('gratrast_move%d%s_%d_%d',rep,Block_Name,channel_no,clust_no)),'fig')
%         saveas(hist_fig,fullfile(apname,sprintf('grathist_move%d%s_%d_%d',rep,Block_Name,channel_no,clust_no)),'fig');
%         %saveas(fft_fig,fullfile(pname,sprintf('gratfft_move%d%s_%d_%d',rep,Block_Name,channel_no,clust_no)),'fig');
        cell_n
        close(rast_fig);
        clear rast_fig;
    end %%% rep
    both_theta=figure
    plot(theta_ind*180/pi,drift_thetatuning(cell_n,:,1)+both_drift_spont(cell_n,1))
    hold on
    plot(theta_ind*180/pi,ones(12,1)*both_drift_spont(cell_n,1),':')
    plot(theta_ind*180/pi,drift_thetatuning(cell_n,:,2)+both_drift_spont(cell_n,2),'g')
    plot(theta_ind*180/pi,ones(12,1)*both_drift_spont(cell_n,2),'g:')
%     saveas(both_theta,fullfile(apname,sprintf('grattheta_move%s_%d_%d',Block_Name,channel_no,clust_no)),'fig')
    
    both_sf=figure
    plot(interp_sfs,drift_sftuning(cell_n,:,1)+both_drift_spont(cell_n,1));
    hold on
    plot(interp_sfs,ones(25,1)*both_drift_spont(cell_n,1),':')
    plot(interp_sfs,ones(25,1)*both_drift_spont(cell_n,2),'g:')
    plot(interp_sfs,drift_sftuning(cell_n,:,2)+both_drift_spont(cell_n,2),'g');
%     saveas(both_sf,fullfile(apname,sprintf('gratsf_move%s_%d_%d',Block_Name,channel_no,clust_no)),'fig')
end
figure
bar(both_drift_spont);

clear event_times_all etimes_old idx_all tm used score c_score channel_times

driftorient_thetawidth
driftorient_A1
driftorient_A2
driftorient_B

osi_est = (driftorient_A1 )./(driftorient_A1 +driftorient_B);

drift_spont
wpref_dog
wbw_dog
peak_value

if use_afile
    
    drift_orients = orients;
    drift_spatfreqs = spatfreqs;
    
    drift_TF = tempfreq;
    
    drift_sponterr = R_err(:,n_cond-1);
    drift_spontvar = R_var(:,n_cond-1);
    
    save(afile_drift, 'drift_orientfreq_all', 'driftorient_wpref', 'driftorient_bw',  'driftorient_theta', 'driftorient_osi','driftorient_dsi','driftorient_A1', 'driftorient_A2', 'driftorient_thetawidth', 'driftorient_B', ...
        'driftorient_F0', 'driftorient_F1', 'driftorient_null', 'drift_spont',  'drift_TF', 'drift_thetatuning','drift_sftuning','drift_sftuning_round', ....
        'drift_peakR','wpref_dog', 'wbw_dog', 'peak_value','-append');
end

invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');