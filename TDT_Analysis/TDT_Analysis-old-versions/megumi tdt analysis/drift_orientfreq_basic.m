function drift_orientfreq_basic
% Matlab code for reading from TTank for drifting or counterphase gratings
% plots histgrams and rasters and performs fourier analysis
%function gratfreq_cluster
% Uses clustering information from cluster_tetrode.m and select_units
% cmn 06-06, based on code by Jianhua Cang 06-27-03
%%% includes two repetitions, for different movement conditions
%%% if no movement, set n_reps =1;

%%%% parameter list
orients = 0:30:330                          %%% orientation list
spatfreqs = [.01 .02 .04 .08 .16 .32];      %%%% spatial frequency list
n_rows=length(orients);
n_col=length(spatfreqs);
stim_duration =1.5;
wait_duration = .5;  %duration after each stimulus (for calculating spontaneous rate)
blank_interval = 0.1; %% length of time after stimulus not to use in calculating spontaneous rate
fft_int = .05;  %%% size of bins to use for fourier analyis (in sec)
tempfreq = 2;    %%% temporal frequency of stimulus
blank_stim = 1; %%% is there an extra stimulus to measure spontaneous
full_field = 1;  %%% is there a full-field flicker?
n_reps = 1;


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

if n_reps>1
    [tsamp vsmooth] = getBlockVelocity(Tank_Name,Block_Name);
    thresh_velocity = 1;
end

Event_Name_Snip='Snip'
Event_Name_Wave='PDec'
Sample_Interval=0.04096 % 24414.0625Hz
Sample_Number_Snip=64
Dec_Factor=32; %
plot_duration = stim_duration + wait_duration;

hist_int = 0.05;
hist_range=[0:hist_int:plot_duration];
axis_range=[0 plot_duration 0 25];
max_events=100000;
ch_map = 1:16;

[afname, apname] = uigetfile('*.mat','analysis data');
if fname~=0
    afile_drift = fullfile(apname,afname);
    save(afile_drift, 'afile_drift','apname','-append'); %%% to make sure it's not overwritten
    load(afile_drift);
    use_afile=1;
else
    use_afile=0;
end
clear drift_orientfreq_all

%%% set time based on first and last timepoints in clustered data
invoke(TTX,'CreateEpocIndexing');
event_times_all = event_times_all-(block-1)*10^5;
times = event_times_all(event_times_all>0 & event_times_all<10^5);

times = times(times>0);
time1 = min(times)
time2= max(times)
MyEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', time1,time2, 5000);
Screen_Plot_y = [axis_range(3):0.1:axis_range(4)]; % to mark the bar appearing and disappearing time

% set number of conditions and display setup (generally rows = orientation, columns = frequency)

n_cond=n_rows*n_col;
if blank_stim
    n_cond=n_cond+1;
end
if full_field
    n_cond=n_cond+1;
end

%printfig = input('print ? ');
printfig=0;
for cell_n = 1:size(cells,1)
    channel_no = cells(cell_n,1)
    clust_no = cells(cell_n,2)
    channel_times =squeeze(event_times_all(channel_no,:));

    times = channel_times(channel_times>0 & channel_times<10^5);
    time1 = min(times);
    time2=max(times);

    for rep =1:n_reps
        hist_fig = figure('Name',sprintf('psth ch:%d cl:%d',channel_no,clust_no));
        rast_fig = figure( 'Name',sprintf('rasters ch:%d cl:%d',channel_no,clust_no));
        spont_full_rast=figure( 'Name',sprintf('spont / flicker ch:%d cl:%d',channel_no,clust_no));
        for cond = [n_cond-2 n_cond-1 0:n_cond-3]
            clear Epocs numtrials Epocs_TS N Spike_TS index TS_xTrg Spike_Timing;
            Epocs=find(MyEpocs(1,:)==cond+1);
            numtrials = length(Epocs);
            Epocs_TS = MyEpocs(2,Epocs);

            invoke(TTX,'ResetFilters');
            ecode= invoke(TTX, 'StringToEvCode', 'xTrg'); % to convert string xTrig to code
            invoke(TTX,'SetFilter', ecode, 69, cond+1, 0);
            % 69 means equal to: 'xTrg=orientation', the last parameter not used for '69'
            N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_Snip, ch_map(channel_no), 0, ...
                max(time1,Epocs_TS(1)),time2,  'FILTERED');
            if (N==max_events)
                warning('max number of events aquired');
            end
            Spike_TS = invoke(TTX, 'ParseEvInfoV', 0, N, 6);  %   6  = Time Stamp

            %%% match up spike times for this condition with spike times from clustering
            %%% to identify cluster index for each spike
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


            % Now for all events, we find out when the xTrig was:
            index = zeros(1,N);
            TS_xTrg=index;
            for i = 1:size(Epocs_TS,2);
                if i<size(Epocs_TS,2)
                    epochSpikes = find(Spike_TS>=Epocs_TS(i) & Spike_TS<Epocs_TS(i+1));
                else
                    epochSpikes = find(Spike_TS>=Epocs_TS(i) & Spike_TS<Epocs_TS(i) + stim_duration+wait_duration+0.05);
                end
                index(epochSpikes)=i;
                TS_xTrg(epochSpikes)=Epocs_TS(i);
            end

            Spike_Timing=Spike_TS-TS_xTrg;
            numtrials = max(index);
            % use if numtrails comes up empty due to no activity on a
            % tetrode
            %numtrials = size(Epocs_TS,2);

            Spike_Timing=Spike_TS-TS_xTrg;

            if n_reps>1
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
                    notenough(cond+1)=1;

                else
                    notenough(cond+1)=0;
                end
                speed(cond+1)=trial_velocity(trials(1));

                for i = 1:numtrials;
                    if isempty(find(trials==i))
                        newtrial(i)=0;
                    else
                        newtrial(i) = find(trials==i);
                    end
                end
                numtrials = sum(usedtrial);

                %%% only keep spikes that are in the desired cluster
                Spike_Timing = Spike_Timing(find((clust_idx==clust_no)&(usedtrial(index))));
                Spike_TS=Spike_TS(find((clust_idx==clust_no)&(usedtrial(index))));
                index=index(find((clust_idx==clust_no)&(usedtrial(index))));
                index = newtrial(index);  %%% get rid of unused trials

                %%%% end movememnt specific
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            else  %%% no movement
                Spike_Timing = Spike_Timing(find((clust_idx==clust_no)));
                Spike_TS=Spike_TS(find((clust_idx==clust_no)));
                index=index(find((clust_idx==clust_no)));
            end

            title_text=['Orientation: ' num2str(cond*45)];
            stim_duration
            numtrials
            R(cell_n,cond+1) = sum(Spike_Timing<=stim_duration)/(stim_duration*numtrials);
            %             if rep==1 & numtrials==1
            %                 R(cell_n,cond+1)=0;
            %             end
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
                    title(title_text,'FontSize',8);
                end
                if cond==1
                    title_text = sprintf(' ch%d c%d',channel_no,clust_no);
                    title(title_text,'FontSize',8);
                end
            end
            %% histograms
            rate_hist = hist(Spike_Timing, hist_range)/(hist_int*numtrials);
            if cond<n_rows*n_col
                figure(hist_fig);
                subplot(n_rows,n_col,cond+1);
                bar(hist_range, rate_hist);  hold on;
                plot([0 stim_duration],[spont_avg spont_avg],'g');
                axis(axis_range);
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
                if cond==0
                    title_text = sprintf('%s',Tank_Name);
                    title(title_text,'FontSize',8);

                elseif cond==1
                    title_text = sprintf('ch%d c%d',channel_no,clust_no);
                    title(title_text,'FontSize',8);
                end
            end

            if cond>=n_rows*n_col & full_field     %%% blank frame or full-field, but only plot if full-field is on
                %
                %                     %%% don't display hist for flicker/spont - too many figures
                %                     figure(spont_full)
                %                     subplot(2,1,cond-n_rows*n_col +1);
                %                     bar(hist_range, rate_hist);  hold on;
                %                     axis(axis_range);
                %                     set(gca,'XTickLabel',[])
                %                     set(gca,'YTickLabel',[])
                %                     title(title_text);

                figure(spont_full_rast);
                if cond == n_cond-2
                    spont_avg = sum(Spike_Timing<=stim_duration)/(stim_duration*numtrials);
                end
                subplot(2,1,cond-n_rows*n_col +1);
                hold on; set(gca, 'yDir','reverse');
                plot (Spike_Timing, index, '.k', 'MarkerSize',8);
                axis([0 plot_duration 0 numtrials+1]);
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
                if cond==n_cond-2
                    title_text = sprintf('spont %s',Tank_Name);
                    title(title_text,'FontSize',8);
                end
                if cond==n_cond-1
                    title_text = sprintf('flicker ch%d c%d',channel_no,clust_no);
                    title(title_text,'FontSize',8);
                end

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

            spont(cell_n,cond+1) = sum((Spike_Timing>(stim_duration+blank_interval))&(Spike_Timing<(stim_duration+wait_duration)))/(numtrials*(wait_duration-blank_interval));

            F0(cell_n,cond+1) = fft_data(1);
            F1(cell_n,cond+1) = 2*fft_data(1+tempfreq/freq_int);  %%% double to count both pos & neg frequency
            F2(cell_n,cond+1) = 2*fft_data(1+2*tempfreq/freq_int);

        end  %% cond

        title_text = sprintf('channel %d cluster %d',channel_no, clust_no);

        if blank_stim
            drift_spont(cell_n,rep)=R(cell_n,n_rows*n_col+1);
        else
            spont_avg = mean(spont(cell_n,:))  %%% if no blank frame, average over inter-stimulus interval of all conditions
        end

        plotcolor = 'bgrcmykbgr';
        tuning_fig = figure('Name', sprintf('All tuning ch:%d cl:%d',channel_no,clust_no));
        %for f= 1:n_col;
        for f = 1:n_col
            plot(orients, R(cell_n,f:n_col:f+n_col*(n_rows-1))-spont_avg,plotcolor(f));
            sftuning_all(:,f,rep) = R(cell_n,f:n_col:f+n_col*(n_rows-1));
            hold on;
        end
        legend('.01 cpd','.02cpd','.04cpd','.08cpd','.16cpd','.32cpd')
        xlabel('deg')
        ylabel('sp / sec')
        %% calculate tuning parameters

        orientfreq = reshape(R(cell_n,1:n_col*n_rows),n_col,n_rows)'-spont_avg;
        orient_tuning_all = mean(orientfreq,2);
        %             figure
        %             plot(orient_tuning_all);
        if max(orient_tuning_all)>abs(min(orient_tuning_all))
            [max_resp pref_orient(cell_n)] = max(orient_tuning_all);
        else
            [max_resp pref_orient(cell_n)] = min(orient_tuning_all);
        end
        freq_tuning(cell_n,:) = orientfreq(pref_orient(cell_n),:);
        %            figure
        %             subplot(2,1,1);
        %             plot(freq_tuning(cell_n,:));
        if max(freq_tuning)>abs(min(freq_tuning))
            [max_resp pref_freq(cell_n)] = max(freq_tuning(cell_n,:));
        else
            [max_resp pref_freq(cell_n)] = min(freq_tuning(cell_n,:));
        end
        orient_tuning(cell_n,:) = orientfreq(:,pref_freq(cell_n))';

        n_orient = n_rows;
        orientfreq(:,2:n_col+1) = reshape(R(cell_n,1:n_col*n_rows),n_col,n_rows)'-spont_avg;
        orientfreq(:,1) = R(cell_n,n_col*n_rows+2)-spont_avg;  %%% full-field (w=0)

        n_w = n_col+1;
        orient_tuning_mn = mean(orientfreq(:,2:n_col+1),2)';

        theta_ind = orients*pi/180;
        [theta_pref OSI A1 A2 w B yfit] = fit_tuningcurve_basic(orient_tuning_mn,theta_ind);

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

        theta_pref_row = round(theta_pref/(orients(2)*pi/180))+1;
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
        if w_pref>n_col
            w_pref=n_col;
        end
        w_pref   %%% this is the index (i.e. column) of preferred w
        w_bw = wbw_dog(cell_n,rep)

        title(title_text);

        %             figure
        %             plot(w_tuning_curve);
        %             title('w tuning curve');
        if max(w_tuning_curve)<abs(min(w_tuning_curve)) & abs(min(w_tuning_curve))>2
            [amp w_pref] = min(w_tuning_curve)
        else
            [amp w_pref] = max(w_tuning_curve)
        end
        w_pref=interp_sfs(w_pref)
        %%% fit orientation tuning curve at optimal SF
        w_pref_col = ones(size(1:n_orient)).*w_pref;
        theta_tuning = interp2(w,theta,orientfreq,w_pref_col,theta_ind);
        sprintf('%d',w_pref)

        [driftorient_theta(cell_n,rep) driftorient_osi(cell_n,rep) driftorient_A1(cell_n,rep) driftorient_A2(cell_n,rep) ...
            driftorient_thetawidth(cell_n,rep) driftorient_B(cell_n,rep) driftorient_null(cell_n,rep) drift_curvefit(cell_n,rep,:)]= fit_tuningcurve_basic(theta_tuning,theta_ind);

        
        [drift_OSI(cell_n,rep) drift_thetawidth(cell_n,rep) drift_peak(cell_n,rep)] = ...
            calculate_tuning(driftorient_A1(cell_n,rep),driftorient_A2(cell_n,rep),driftorient_B(cell_n,rep),driftorient_thetawidth(cell_n,rep));

        
        i=sqrt(-1);
        mu = (sum(theta_tuning.*exp(i*theta_ind)))/sum(theta_tuning);
        if isnan(mu);
            mu=0;
        end

        driftorient_dsi(cell_n) = abs(mu);

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

        title_text = sprintf('%s ch%d cl%d',Tank_Name, channel_no,clust_no);

        if printfig
            print(hist_fig);
        end


        % saveas(tuning_fig,fullfile(pname,sprintf('grattuning_move%d%s_%d_%d',rep,Block_Name,channel_no,clust_no)),'fig')
        saveas(rast_fig,fullfile(apname,sprintf('gratrast_move%d%s_%d_%d',rep,Block_Name,channel_no,clust_no)),'fig')
        saveas(hist_fig,fullfile(apname,sprintf('grathist_move%d%s_%d_%d',rep,Block_Name,channel_no,clust_no)),'fig');
        %saveas(fft_fig,fullfile(pname,sprintf('gratfft_move%d%s_%d_%d',rep,Block_Name,channel_no,clust_no)),'fig');

    end %%% rep
    both_theta=figure('Name',sprintf('orient tuning ch:%d cl:%d',channel_no,clust_no))
    plot(theta_ind*180/pi,drift_thetatuning(cell_n,:,1)+drift_spont(cell_n,1),'b')
    hold on
    plot(theta_ind*180/pi,ones(n_rows,1)*drift_spont(cell_n,1),'b:')
    plot(theta_ind*180/pi,squeeze(drift_curvefit(cell_n,1,1:length(orients)))+drift_spont(cell_n,1),'g:')
   legend('data','spont','fit')
    if n_reps>1
        plot(theta_ind*180/pi,drift_thetatuning(cell_n,:,2)+drift_spont(cell_n,2),'r')
        plot(theta_ind*180/pi,ones(n_rows,1)*drift_spont(cell_n,2),'r:')
        plot(theta_ind*180/pi,drift_curvefit(cell_n,2,:)+drift_spont(cell_n,1),'b:')
    end
xlabel('orientation (deg)')
ylabel('sp / sec')
    saveas(both_theta,fullfile(apname,sprintf('grattheta_move%s_%d_%d',Block_Name,channel_no,clust_no)),'fig')

    %         both_sf=figure
    %         plot(interp_sfs,drift_sftuning(cell_n,:,1)+drift_spont(cell_n,1));
    %         hold on
    %         plot(interp_sfs,ones(25,1)*drift_spont(cell_n,1),':')
    %         plot(interp_sfs,ones(25,1)*drift_spont(cell_n,2),'g:')
    %         plot(interp_sfs,drift_sftuning(cell_n,:,2)+drift_spont(cell_n,2),'g');
    %         saveas(both_sf,fullfile(apname,sprintf('gratsf_move%s_%d_%d',Block_Name,channel_no,clust_no)),'fig')
end

clear event_times_all etimes_old idx_all tm used score c_score channel_times


if use_afile
    drift_orients = orients;
    drift_spatfreqs = spatfreqs;
    drift_TF = tempfreq;
    
        save(afile_drift, 'drift_orientfreq_all', 'driftorient_wpref', 'driftorient_bw',  'driftorient_theta', 'driftorient_osi','driftorient_dsi','driftorient_A1', 'driftorient_A2', 'driftorient_thetawidth', 'driftorient_B', ...
            'driftorient_F0', 'driftorient_F1', 'driftorient_null', 'drift_spont',  'drift_TF', 'drift_thetatuning','drift_sftuning','drift_sftuning_round', ....
            'wpref_dog', 'wbw_dog', 'drift_OSI','drift_thetawidth','drift_peak','-append');
end

drift_spont
drift_peak
drift_OSI
drift_thetawidth

invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');