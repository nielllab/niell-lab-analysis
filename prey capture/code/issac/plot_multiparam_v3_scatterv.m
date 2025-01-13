% approach freq, intercept given approach, capture given intercept

addpath('f:\prey_capture_lkhd_analysis\code')  
addpath('f:\code')  
addpath('f:\code\hline_vline\')
addpath('f:\code\povilaskarvelis\daboxplot\')
addpath('f:\code\densityScatterChart-1.2.0.0\')

% load 'app_d'
% load('f:\prey_capture_lkhd_analysis\d_save\app_d_v9.mat')
% load('f:\prey_capture_lkhd_analysis\d_save\app_d_v11.mat')
load('f:\prey_capture_lkhd_analysis\d_save\app_d_v15.mat')
% load corrected capture time 'df_meta'
load('f:\prey_capture_lkhd_analysis\d_save\df_meta_v1.mat')

% condition list
cond_list = {'Wno','Lno','Hno','Wsw','Wlb','Wsb','Lsb','Hsb'};

% create output var
app_len = []; % lengths of each approach
app_freq = []; % frequency of approach, #/time
app_num = []; % number of approaches
app_perc = []; % approach as percentage of total time
app_tot = []; % total time of approach
int_perc = []; % intercept as percentage of total time
int_num = []; % number of intercepts
capt_bin = []; % capture binary
capt_time = []; % capture time (fr)
success_app = []; % success rate per app (success/attempt #)
success_int = []; % success rate per int (success/attempt #)
int_per_app = []; % probability of intercept given approach
int_per_app_stat = []; % grab the statistic for int per app
for icond = 1:8
    for ilaser = 1:2
        app_len{icond,ilaser} = [];
        app_freq{icond,ilaser} = [];
        app_num{icond,ilaser} = [];
        app_tot{icond,ilaser} = [];
        app_perc{icond,ilaser} = [];
        int_perc{icond,ilaser} = [];
        int_num{icond,ilaser} = [];
        capt_bin{icond,ilaser} = [];
        capt_time{icond,ilaser} = [];
        success_app{icond,ilaser} = [];
        success_int{icond,ilaser} = [];
        int_per_app{icond,ilaser} = [];
        int_per_app_stat{icond,ilaser} = [];
    end
end

% icond = 1;
% targetcond = cond_list{icond};
% targetlaser = 'n';
targetexpctrl = 'Exp';

% first, histograms of lengths of approach distributions
for ifile = 1:size(app_h,2)

    % note that for 'app_d_v8.mat', app_h has the correct index, everything
    % else is +1
    fname = app_h{ifile}.fname;
    fdate = app_h{ifile}.fdate; % only used to locate target_row in df_meta
    fcond = app_h{ifile}.fcond;
    flaser = app_h{ifile}.flaser;
    ftrial = app_h{ifile}.ftrial;

    % if strcmp(fcond,targetcond) == 1 && strcmp(flaser,targetlaser) == 1 

    cond_loc = [];
    for icond = 1:8
        if strcmp(fcond,cond_list{icond})
            cond_loc = icond;
        end
    end
    if isempty(cond_loc)
        disp(['Could not find the correct condition for this trial, something is wrong: ifile #' num2str(ifile)])
    end

    laser_loc = [];
    if strcmp(flaser,'n')
        laser_loc = 1;
    elseif strcmp(flaser,'y')
        laser_loc = 2;
    end

    % get capture time (if any)
    capt_binary = [];
    capt_t = [];
    % find df_meta correspondence
    df_loc = [];
    for idf = 2:size(df_meta,1)
        df_name = df_meta{idf,17};
        df_date = num2str(df_meta{idf,18});
        if length(df_date) == 5
            df_date = ['0' df_date];
        end
        df_cond = df_meta{idf,9};
        df_trial = num2str(df_meta{idf,16});
        df_laser = df_meta{idf,13};
        if strcmp(df_laser,'True')
            df_laser = 'y';
        elseif strcmp(df_laser,'False')
            df_laser = 'n';
        else
            disp('What else is there here other than laser-on vs off???')
        end

        if strcmp(df_name,fname) && strcmp(df_date,fdate) && strcmp(df_cond,fcond) && strcmp(df_trial,ftrial) && strcmp(df_laser,flaser)
            df_loc = idf;
            expctrl = df_meta{df_loc,10};
        end
    end

    if strcmp(expctrl,targetexpctrl)


        disp(['Currently: ' num2str(ifile) ' of ' num2str(size(app_h,2)) ' ... ' fname ', ' fdate ', ' fcond ', ' flaser ', ' ftrial])

        capt_t = df_meta{df_loc,2}; capt_t = round(capt_t*60); % in unit of frames
        if isempty(capt_t) == 1
            capt_binary = 0;
        elseif isempty(capt_t) == 0
            if capt_t < 30*60 % it's in frames now
                capt_binary = 1;
            else
                capt_binary = 0;
            end
        end
    
        % total length of the experiment; the capture cut-off is not previously
        % built in (doing it here)
        % use corrected capture time 'df_meta' .mat file
        % add 60 frames (1s) when capture
        file_id = ifile + 0; % correct for 'app_d_v8' here by +1
    
        if capt_binary == 0 % no capture
            trial_len = 30*60; % frames
        elseif capt_binary == 1
            trial_len = capt_t;
            param_len = size(app_p{file_id},2);
            if (trial_len+60) < param_len
                trial_len = trial_len + 60;
            elseif (trial_len+60) > param_len
                trial_len = param_len;
            end
        end
    
    
        % calculate the number of approaches and the lengths of each approaches
        tmp_appd = app_d{file_id};
        if isempty(tmp_appd) == 0
            if capt_binary == 0 % no capture
                thrs = 30*60;
                exclude_list = find(tmp_appd(:,1) > thrs);
                if isempty(exclude_list) == 0
                    tmp_appd(exclude_list,:) = [];
                    if isempty(tmp_appd) == 0
                        if tmp_appd(end,2) > thrs
                            tmp_appd(end,2) = thrs;
                        end
                    end 
                end
            elseif capt_binary == 1
                param_len = size(app_p{file_id},2);
                thrs = capt_t + 60; % 60 is a cushion
                % make sure the cushion + capt_t is not longer than the actual data
                if thrs > param_len
                    thrs = param_len;
                end
                exclude_list = find(tmp_appd(:,1) > thrs);
                if isempty(exclude_list) == 0
                    tmp_appd(exclude_list,:) = [];
                    if isempty(tmp_appd) == 0
                        if tmp_appd(end,2) > thrs
                            tmp_appd(end,2) = thrs;
                        end
                    end
                end
            end
    
            tmp_len = tmp_appd(:,2) - tmp_appd(:,1) + 1; % +1 is to include the subtracted frame
            num_app = length(tmp_len); % # of app
            tot_app = sum(tmp_len); % # (frames) consumed by app
            avg_app = mean(tmp_len); % average duration of app
        
            p_spd = app_p{file_id}(1,:);
            p_ang = app_p{file_id}(2,:);
            p_d2c = app_p{file_id}(3,:);
            app_int = [];
            for iapp = 1:size(tmp_appd,1)
                st = tmp_appd(iapp,1);
                ed = tmp_appd(iapp,2);
                ed1 = ed;
                if ed == size(app_p{file_id},2)
                    ed;
                else
                    ed = ed + 1;
                end
        
                tmp_d2c = p_d2c(st:ed);
                f4  = find(tmp_d2c < 4);
                if isempty(f4) == 0
                    app_int = [app_int; 1]; % intercept happened
                elseif isempty(f4) == 1
                    app_int = [app_int; 0]; % intercept did not happen
                end
            end
            
            prob_int_given_app = sum(app_int)/length(app_int);
        
            num_int = sum(app_int); % # of int
        
            % find the approach that eventually led to the capture
            if capt_binary == 1
                app_list = find(tmp_appd(:,1) < capt_t);
                app_loc = app_list(end);
                st = tmp_appd(app_loc,1);
                ed = tmp_appd(app_loc,2);
        
                app_c_tot = ed-st+1; % length
        
                % success rate (per app)
                app_success = 1/size(tmp_appd,1);
                int_success = 1/num_int;
            else
                app_success = 0;
                int_success = 0;
            end
        
            perc_app = tot_app/trial_len;
            tot_int_dur = 0;
            for iint = 1:length(app_int)
                if app_int(iint) == 1
                    st = tmp_appd(iint,1);
                    ed = tmp_appd(iint,2);
                    tot_int_dur = tot_int_dur + ed - st + 1;
                end
            end
            perc_int = tot_int_dur/trial_len;
            freq_app = num_app/trial_len;
            freq_int = num_int/trial_len;
        
            % save out
            for icond = cond_loc
                for ilaser = laser_loc
                    app_len{icond,ilaser} = [app_len{icond,ilaser}; avg_app];
                    app_tot{icond,ilaser} = [app_tot{icond,ilaser}; tot_app];
                    app_num{icond,ilaser} = [app_num{icond,ilaser}; num_app];
                    app_perc{icond,ilaser} = [app_perc{icond,ilaser}; perc_app];
                    app_freq{icond,ilaser} = [app_freq{icond,ilaser}; freq_app];
        
                    int_perc{icond,ilaser} = [int_perc{icond,ilaser}; perc_int];
                    int_num{icond,ilaser} = [int_num{icond,ilaser}; num_int];
                    if isempty(capt_t) == 1
                        capt_time{icond,ilaser} = [capt_time{icond,ilaser}; NaN];
                        capt_bin{icond,ilaser} = [capt_bin{icond,ilaser}; capt_binary];
                    elseif isempty(capt_t) == 0
                        capt_time{icond,ilaser} = [capt_time{icond,ilaser}; capt_t];
                        capt_bin{icond,ilaser} = [capt_bin{icond,ilaser}; capt_binary]; % capt_binary keeps track whether the capture, even if it happened, happened prior to 30s (which is 1, after is still 0)
                    end
        
                    success_app{icond,ilaser} = [success_app{icond,ilaser}; app_success];
                    success_int{icond,ilaser} = [success_int{icond,ilaser}; int_success];

                    int_per_app{icond,ilaser} = [int_per_app{icond,ilaser}; prob_int_given_app];
                    if isempty(capt_t) == 1
                        int_per_app_stat{icond,ilaser} = [int_per_app_stat{icond,ilaser}; num_app num_int 1];
                    elseif isempty(capt_t) == 0
                        int_per_app_stat{icond,ilaser} = [int_per_app_stat{icond,ilaser}; num_app num_int 0];
                    end
                end
            end
    
        elseif isempty(tmp_appd) == 1
    
            % save out
            for icond = cond_loc
                for ilaser = laser_loc
                    app_len{icond,ilaser} = [app_len{icond,ilaser}; 0];
                    app_tot{icond,ilaser} = [app_tot{icond,ilaser}; 0];
                    app_num{icond,ilaser} = [app_num{icond,ilaser}; 0];
                    app_perc{icond,ilaser} = [app_perc{icond,ilaser}; 0];
                    app_freq{icond,ilaser} = [app_freq{icond,ilaser}; 0];
        
                    int_perc{icond,ilaser} = [int_perc{icond,ilaser}; 0];
                    int_num{icond,ilaser} = [int_num{icond,ilaser}; 0];
                    if isempty(capt_t) == 1
                        capt_time{icond,ilaser} = [capt_time{icond,ilaser}; NaN];
                        capt_bin{icond,ilaser} = [capt_bin{icond,ilaser}; 0];
                    elseif isempty(capt_t) == 0
                        capt_time{icond,ilaser} = [capt_time{icond,ilaser}; capt_t];
                        capt_bin{icond,ilaser} = [capt_bin{icond,ilaser}; capt_binary]; % capt_binary keeps track whether the capture, even if it happened, happened prior to 30s (which is 1, after is still 0)
                    end
        
                    success_app{icond,ilaser} = [success_app{icond,ilaser}; NaN];
                    success_int{icond,ilaser} = [success_int{icond,ilaser}; NaN];

                    int_per_app{icond,ilaser} = [int_per_app{icond,ilaser}; NaN];
                    if isempty(capt_t) == 1
                        int_per_app_stat{icond,ilaser} = [int_per_app_stat{icond,ilaser}; 0 0 1];
                    elseif isempty(capt_t) == 0
                        int_per_app_stat{icond,ilaser} = [int_per_app_stat{icond,ilaser}; 0 0 0];
                    end
                end
            end
        end

    end
    % end
end


    
%% Figures

% cd('F:\prey_capture_lkhd_analysis\figure\2024_11_18\') % save here (save folder sav fldr)
cd('F:\prey_capture_lkhd_analysis\figure\2025_01_06c\')

% frequency of approach (count as absolute # divided by time)
f = figure(1829); clf
f.Position = [114 227 1795 625];
h = []; xlimcmp = []; ylimcmp = [];
rho = []; pval = [];
pev = [];
for icond = 1:8
    for ilaser = 1:2
        pev{icond,ilaser} = [];
    end
end
for icond = 1:8
    sav_faildat = [];
    for ilaser = 1:2

        sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
        if ilaser == 1
            cla
        end

        xp1 = app_freq{icond,ilaser}/60;
        yp1 = capt_time{icond,ilaser}/60;
        yp1_nan = isnan(yp1);
        yp1_nanone = find(yp1_nan == 1);
        yp1_over = find(yp1>(1800/60));
        yp1_clist = [yp1_nanone; yp1_over];

        xp1_success = xp1; xp1_success(yp1_clist) = [];
        yp1_success = yp1; yp1_success(yp1_clist) = [];

        xp1_fail = xp1(yp1_clist);
        yp1_fail = yp1(yp1_clist);
        sav_faildat{ilaser} = xp1_fail;

        % regression:
        pev{icond,ilaser} = polyfit(yp1_success,xp1_success,1);
        [rho(icond,ilaser),pval(icond,ilaser)] = corr(yp1_success,xp1_success);

        h{icond,ilaser} = plot(yp1_success,xp1_success,'o'); hold on
        if ilaser == 1
            h{icond,ilaser}.Color = [0 0 0];
            h{icond,ilaser}.MarkerFaceColor = [0 0 0];
            h{icond,ilaser}.MarkerSize = 3;
        elseif ilaser == 2
            h{icond,ilaser}.Color = [0 0.75 0.5];
            h{icond,ilaser}.MarkerFaceColor = [0 0.75 0.5];
            h{icond,ilaser}.MarkerSize = 3;
        end
        
        ax1 = gca;
        xlimcmp = [xlimcmp; ax1.XLim];
        ylimcmp = [ylimcmp; ax1.YLim];

        sub2 = subplot(2,12,icond*3);
        if ilaser == 1
            cla
        end
        if ilaser == 2

            hbox = daboxplot(sav_faildat,'scatter',2,'whiskers',0,'boxalpha',0.7);
            ax2 = gca;
            ax2.XTickLabel = {'Loff','Lon'};

            hbox.bx(1).FaceColor = [0 0 0];
            hbox.bx(1).FaceAlpha = 0.3;

            hbox.bx(2).FaceColor = [0 0.75 0.5];
            hbox.bx(2).FaceAlpha = 0.3;
            hbox.sc(2).MarkerEdgeColor = [1 1 1];
            hbox.sc(2).MarkerFaceColor = [0 0.75 0.5];
        end
        

    end
end
minxlim = min(xlimcmp(:,1));
maxxlim = max(xlimcmp(:,2));
minylim = min(ylimcmp(:,1));
maxylim = max(ylimcmp(:,2));
plotx = minxlim:(maxxlim-minxlim)/100:maxxlim;
for icond = 1:8
    sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
    xlim([minxlim maxxlim])
    ylim([minylim maxylim])

    for ilaser = 1:2
        vala = pev{icond,ilaser}(1);
        valint = pev{icond,ilaser}(2);
        ploty = plotx*vala + valint;
        ptmp = plot(plotx,ploty); hold on
        if ilaser == 1
            ptmp.Color = [0 0 0];
            % ptmp.LineWidth = 2;
            ptmp.LineStyle = ':'
        elseif ilaser == 2
            ptmp.Color = [0 0.75 0.5];
            ptmp.LineStyle = '--';
        end
    end

    title({cond_list{icond},['r=' num2str(rho(icond,1)) ', p=' num2str(pval(icond,1))],['r=' num2str(rho(icond,2)) ', p=' num2str(pval(icond,2))]})


    if icond == 5 || icond == 6 || icond == 7 || icond == 8
        xlabel('capt time (s)')
    end

    if icond == 1 || icond == 5
        ylabel('app freq (#/time(s))')
    end


    sub2 = subplot(2,12,icond*3);
    ylim([minylim maxylim])
end

F = getframe(f);
imwrite(F.cdata, 'scatter_app_freq.png', 'png')
close








% percentage of time spent in approach
f = figure(1849); clf
f.Position = [114 227 1795 625];
h = []; xlimcmp = []; ylimcmp = [];
rho = []; pval = [];
pev = [];
for icond = 1:8
    for ilaser = 1:2
        pev{icond,ilaser} = [];
    end
end
for icond = 1:8
    sav_faildat = [];
    for ilaser = 1:2

        sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
        if ilaser == 1
            cla
        end

        xp1 = app_perc{icond,ilaser};
        yp1 = capt_time{icond,ilaser}/60;
        yp1_nan = isnan(yp1);
        yp1_nanone = find(yp1_nan == 1);
        yp1_over = find(yp1>(1800/60));
        yp1_clist = [yp1_nanone; yp1_over];

        xp1_success = xp1; xp1_success(yp1_clist) = [];
        yp1_success = yp1; yp1_success(yp1_clist) = [];

        xp1_fail = xp1(yp1_clist);
        yp1_fail = yp1(yp1_clist);
        sav_faildat{ilaser} = xp1_fail;

        % regression:
        pev{icond,ilaser} = polyfit(yp1_success,xp1_success,1);
        [rho(icond,ilaser),pval(icond,ilaser)] = corr(yp1_success,xp1_success);

        h{icond,ilaser} = plot(yp1_success,xp1_success,'o'); hold on
        if ilaser == 1
            h{icond,ilaser}.Color = [0 0 0];
            h{icond,ilaser}.MarkerFaceColor = [0 0 0];
            h{icond,ilaser}.MarkerSize = 3;
        elseif ilaser == 2
            h{icond,ilaser}.Color = [0 0.75 0.5];
            h{icond,ilaser}.MarkerFaceColor = [0 0.75 0.5];
            h{icond,ilaser}.MarkerSize = 3;
        end
        
        ax1 = gca;
        xlimcmp = [xlimcmp; ax1.XLim];
        ylimcmp = [ylimcmp; ax1.YLim];

        sub2 = subplot(2,12,icond*3);
        if ilaser == 1
            cla
        end
        if ilaser == 2

            hbox = daboxplot(sav_faildat,'scatter',2,'whiskers',0,'boxalpha',0.7);
            ax2 = gca;
            ax2.XTickLabel = {'Loff','Lon'};

            hbox.bx(1).FaceColor = [0 0 0];
            hbox.bx(1).FaceAlpha = 0.3;

            hbox.bx(2).FaceColor = [0 0.75 0.5];
            hbox.bx(2).FaceAlpha = 0.3;
            hbox.sc(2).MarkerEdgeColor = [1 1 1];
            hbox.sc(2).MarkerFaceColor = [0 0.75 0.5];
        end
        

    end
end
minxlim = min(xlimcmp(:,1));
maxxlim = max(xlimcmp(:,2));
minylim = min(ylimcmp(:,1));
maxylim = max(ylimcmp(:,2));
plotx = minxlim:(maxxlim-minxlim)/100:maxxlim;
for icond = 1:8
    sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
    xlim([minxlim maxxlim])
    ylim([minylim maxylim])

    for ilaser = 1:2
        vala = pev{icond,ilaser}(1);
        valint = pev{icond,ilaser}(2);
        ploty = plotx*vala + valint;
        ptmp = plot(plotx,ploty); hold on
        if ilaser == 1
            ptmp.Color = [0 0 0];
            % ptmp.LineWidth = 2;
            ptmp.LineStyle = ':'
        elseif ilaser == 2
            ptmp.Color = [0 0.75 0.5];
            ptmp.LineStyle = '--';
        end
    end

    title({cond_list{icond},['r=' num2str(rho(icond,1)) ', p=' num2str(pval(icond,1))],['r=' num2str(rho(icond,2)) ', p=' num2str(pval(icond,2))]})


    if icond == 5 || icond == 6 || icond == 7 || icond == 8
        xlabel('capt time (s)')
    end

    if icond == 1 || icond == 5
        ylabel('app perc (% of time)')
    end


    sub2 = subplot(2,12,icond*3);
    ylim([minylim maxylim])
end
F = getframe(f);
imwrite(F.cdata, 'scatter_app_perc.png', 'png')
close





% probability of intercept given approach
f = figure(1949); clf
f.Position = [114 227 1795 625];
h = []; xlimcmp = []; ylimcmp = [];
rho = []; pval = [];
pev = [];
for icond = 1:8
    for ilaser = 1:2
        pev{icond,ilaser} = [];
    end
end
for icond = 1:8
    sav_faildat = [];
    for ilaser = 1:2

        sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
        if ilaser == 1
            cla
        end

        xp1 = int_per_app{icond,ilaser};
        yp1 = capt_time{icond,ilaser}/60;
        yp1_nan = isnan(yp1);
        yp1_nanone = find(yp1_nan == 1);
        yp1_over = find(yp1>(1800/60));
        yp1_clist = [yp1_nanone; yp1_over];

        xp1_success = xp1; xp1_success(yp1_clist) = [];
        yp1_success = yp1; yp1_success(yp1_clist) = [];

        xp1_fail = xp1(yp1_clist);
        yp1_fail = yp1(yp1_clist);
        sav_faildat{ilaser} = xp1_fail;

        % regression:
        pev{icond,ilaser} = polyfit(yp1_success,xp1_success,1);
        [rho(icond,ilaser),pval(icond,ilaser)] = corr(yp1_success,xp1_success);

        h{icond,ilaser} = plot(yp1_success,xp1_success,'o'); hold on
        if ilaser == 1
            h{icond,ilaser}.Color = [0 0 0];
            h{icond,ilaser}.MarkerFaceColor = [0 0 0];
            h{icond,ilaser}.MarkerSize = 3;
        elseif ilaser == 2
            h{icond,ilaser}.Color = [0 0.75 0.5];
            h{icond,ilaser}.MarkerFaceColor = [0 0.75 0.5];
            h{icond,ilaser}.MarkerSize = 3;
        end
        
        ax1 = gca;
        xlimcmp = [xlimcmp; ax1.XLim];
        ylimcmp = [ylimcmp; ax1.YLim];

        sub2 = subplot(2,12,icond*3);
        if ilaser == 1
            cla
        end
        if ilaser == 2

            hbox = daboxplot(sav_faildat,'scatter',2,'whiskers',0,'boxalpha',0.7);
            ax2 = gca;
            ax2.XTickLabel = {'Loff','Lon'};

            hbox.bx(1).FaceColor = [0 0 0];
            hbox.bx(1).FaceAlpha = 0.3;

            hbox.bx(2).FaceColor = [0 0.75 0.5];
            hbox.bx(2).FaceAlpha = 0.3;
            hbox.sc(2).MarkerEdgeColor = [1 1 1];
            hbox.sc(2).MarkerFaceColor = [0 0.75 0.5];
        end
        

    end
end
minxlim = min(xlimcmp(:,1));
maxxlim = max(xlimcmp(:,2));
minylim = min(ylimcmp(:,1));
maxylim = max(ylimcmp(:,2));
plotx = minxlim:(maxxlim-minxlim)/100:maxxlim;
for icond = 1:8
    sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
    xlim([minxlim maxxlim])
    ylim([minylim maxylim])

    for ilaser = 1:2
        vala = pev{icond,ilaser}(1);
        valint = pev{icond,ilaser}(2);
        ploty = plotx*vala + valint;
        ptmp = plot(plotx,ploty); hold on
        if ilaser == 1
            ptmp.Color = [0 0 0];
            % ptmp.LineWidth = 2;
            ptmp.LineStyle = ':'
        elseif ilaser == 2
            ptmp.Color = [0 0.75 0.5];
            ptmp.LineStyle = '--';
        end
    end

    title({cond_list{icond},['r=' num2str(rho(icond,1)) ', p=' num2str(pval(icond,1))],['r=' num2str(rho(icond,2)) ', p=' num2str(pval(icond,2))]})


    if icond == 5 || icond == 6 || icond == 7 || icond == 8
        xlabel('capt time (s)')
    end

    if icond == 1 || icond == 5
        ylabel('int per app')
    end


    sub2 = subplot(2,12,icond*3);
    ylim([minylim maxylim])
end
F = getframe(f);
imwrite(F.cdata, 'scatter_int_per_app.png', 'png')
close




% histogram of num of approach
f = figure(1221); clf
f.Position = [114 227 1795 625];
h = []; xlimcmp = []; ylimcmp = [];
rho = []; pval = [];
pev = [];
for icond = 1:8
    for ilaser = 1:2
        pev{icond,ilaser} = [];
    end
end
for icond = 1:8
    sav_faildat = [];
    for ilaser = 1:2

        sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
        if ilaser == 1
            cla
        end

        xp1 = app_num{icond,ilaser};
        yp1 = capt_time{icond,ilaser}/60;
        yp1_nan = isnan(yp1);
        yp1_nanone = find(yp1_nan == 1);
        yp1_over = find(yp1>(1800/60));
        yp1_clist = [yp1_nanone; yp1_over];

        xp1_success = xp1; xp1_success(yp1_clist) = [];
        yp1_success = yp1; yp1_success(yp1_clist) = [];

        xp1_fail = xp1(yp1_clist);
        yp1_fail = yp1(yp1_clist);
        sav_faildat{ilaser} = xp1_fail;

        % regression:
        pev{icond,ilaser} = polyfit(yp1_success,xp1_success,1);
        [rho(icond,ilaser),pval(icond,ilaser)] = corr(yp1_success,xp1_success);

        h{icond,ilaser} = plot(yp1_success,xp1_success,'o'); hold on
        if ilaser == 1
            h{icond,ilaser}.Color = [0 0 0];
            h{icond,ilaser}.MarkerFaceColor = [0 0 0];
            h{icond,ilaser}.MarkerSize = 3;
        elseif ilaser == 2
            h{icond,ilaser}.Color = [0 0.75 0.5];
            h{icond,ilaser}.MarkerFaceColor = [0 0.75 0.5];
            h{icond,ilaser}.MarkerSize = 3;
        end
        
        ax1 = gca;
        xlimcmp = [xlimcmp; ax1.XLim];
        ylimcmp = [ylimcmp; ax1.YLim];

        sub2 = subplot(2,12,icond*3);
        if ilaser == 1
            cla
        end
        if ilaser == 2

            hbox = daboxplot(sav_faildat,'scatter',2,'whiskers',0,'boxalpha',0.7);
            ax2 = gca;
            ax2.XTickLabel = {'Loff','Lon'};

            hbox.bx(1).FaceColor = [0 0 0];
            hbox.bx(1).FaceAlpha = 0.3;

            hbox.bx(2).FaceColor = [0 0.75 0.5];
            hbox.bx(2).FaceAlpha = 0.3;
            hbox.sc(2).MarkerEdgeColor = [1 1 1];
            hbox.sc(2).MarkerFaceColor = [0 0.75 0.5];
        end
        

    end
end
minxlim = min(xlimcmp(:,1));
maxxlim = max(xlimcmp(:,2));
minylim = min(ylimcmp(:,1));
maxylim = max(ylimcmp(:,2));
plotx = minxlim:(maxxlim-minxlim)/100:maxxlim;
for icond = 1:8
    sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
    xlim([minxlim maxxlim])
    ylim([minylim maxylim])

    for ilaser = 1:2
        vala = pev{icond,ilaser}(1);
        valint = pev{icond,ilaser}(2);
        ploty = plotx*vala + valint;
        ptmp = plot(plotx,ploty); hold on
        if ilaser == 1
            ptmp.Color = [0 0 0];
            % ptmp.LineWidth = 2;
            ptmp.LineStyle = ':'
        elseif ilaser == 2
            ptmp.Color = [0 0.75 0.5];
            ptmp.LineStyle = '--';
        end
    end

    title({cond_list{icond},['r=' num2str(rho(icond,1)) ', p=' num2str(pval(icond,1))],['r=' num2str(rho(icond,2)) ', p=' num2str(pval(icond,2))]})


    if icond == 5 || icond == 6 || icond == 7 || icond == 8
        xlabel('capt time (s)')
    end

    if icond == 1 || icond == 5
        ylabel('app num (#)')
    end


    sub2 = subplot(2,12,icond*3);
    ylim([minylim maxylim])
end
F = getframe(f);
imwrite(F.cdata, 'scatter_app_num.png', 'png')
close






% avg approach lengths
f = figure(7811); clf
f.Position = [114 227 1795 625];
h = []; xlimcmp = []; ylimcmp = [];
rho = []; pval = [];
pev = [];
for icond = 1:8
    for ilaser = 1:2
        pev{icond,ilaser} = [];
    end
end
for icond = 1:8
    sav_faildat = [];
    for ilaser = 1:2

        sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
        if ilaser == 1
            cla
        end

        xp1 = app_len{icond,ilaser}/60;
        yp1 = capt_time{icond,ilaser}/60;
        yp1_nan = isnan(yp1);
        yp1_nanone = find(yp1_nan == 1);
        yp1_over = find(yp1>(1800/60));
        yp1_clist = [yp1_nanone; yp1_over];

        xp1_success = xp1; xp1_success(yp1_clist) = [];
        yp1_success = yp1; yp1_success(yp1_clist) = [];

        xp1_fail = xp1(yp1_clist);
        yp1_fail = yp1(yp1_clist);
        sav_faildat{ilaser} = xp1_fail;

        % regression:
        pev{icond,ilaser} = polyfit(yp1_success,xp1_success,1);
        [rho(icond,ilaser),pval(icond,ilaser)] = corr(yp1_success,xp1_success);

        h{icond,ilaser} = plot(yp1_success,xp1_success,'o'); hold on
        if ilaser == 1
            h{icond,ilaser}.Color = [0 0 0];
            h{icond,ilaser}.MarkerFaceColor = [0 0 0];
            h{icond,ilaser}.MarkerSize = 3;
        elseif ilaser == 2
            h{icond,ilaser}.Color = [0 0.75 0.5];
            h{icond,ilaser}.MarkerFaceColor = [0 0.75 0.5];
            h{icond,ilaser}.MarkerSize = 3;
        end
        
        ax1 = gca;
        xlimcmp = [xlimcmp; ax1.XLim];
        ylimcmp = [ylimcmp; ax1.YLim];

        sub2 = subplot(2,12,icond*3);
        if ilaser == 1
            cla
        end
        if ilaser == 2

            hbox = daboxplot(sav_faildat,'scatter',2,'whiskers',0,'boxalpha',0.7);
            ax2 = gca;
            ax2.XTickLabel = {'Loff','Lon'};

            hbox.bx(1).FaceColor = [0 0 0];
            hbox.bx(1).FaceAlpha = 0.3;

            hbox.bx(2).FaceColor = [0 0.75 0.5];
            hbox.bx(2).FaceAlpha = 0.3;
            hbox.sc(2).MarkerEdgeColor = [1 1 1];
            hbox.sc(2).MarkerFaceColor = [0 0.75 0.5];
        end
        

    end
end
minxlim = min(xlimcmp(:,1));
maxxlim = max(xlimcmp(:,2));
minylim = min(ylimcmp(:,1));
maxylim = max(ylimcmp(:,2));
plotx = minxlim:(maxxlim-minxlim)/100:maxxlim;
for icond = 1:8
    sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
    xlim([minxlim maxxlim])
    ylim([minylim maxylim])

    for ilaser = 1:2
        vala = pev{icond,ilaser}(1);
        valint = pev{icond,ilaser}(2);
        ploty = plotx*vala + valint;
        ptmp = plot(plotx,ploty); hold on
        if ilaser == 1
            ptmp.Color = [0 0 0];
            % ptmp.LineWidth = 2;
            ptmp.LineStyle = ':'
        elseif ilaser == 2
            ptmp.Color = [0 0.75 0.5];
            ptmp.LineStyle = '--';
        end
    end

    title({cond_list{icond},['r=' num2str(rho(icond,1)) ', p=' num2str(pval(icond,1))],['r=' num2str(rho(icond,2)) ', p=' num2str(pval(icond,2))]})


    if icond == 5 || icond == 6 || icond == 7 || icond == 8
        xlabel('capt time (s)')
    end

    if icond == 1 || icond == 5
        ylabel('app len (time(s))')
    end


    sub2 = subplot(2,12,icond*3);
    ylim([minylim maxylim])
end
F = getframe(f);
imwrite(F.cdata, 'scatter_app_len.png', 'png')
close









% total app time
f = figure(5648); clf
f.Position = [114 227 1795 625];
h = []; xlimcmp = []; ylimcmp = [];
rho = []; pval = [];
pev = [];
for icond = 1:8
    for ilaser = 1:2
        pev{icond,ilaser} = [];
    end
end
for icond = 1:8
    sav_faildat = [];
    for ilaser = 1:2

        sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
        if ilaser == 1
            cla
        end

        xp1 = app_tot{icond,ilaser}/60;
        yp1 = capt_time{icond,ilaser}/60;
        yp1_nan = isnan(yp1);
        yp1_nanone = find(yp1_nan == 1);
        yp1_over = find(yp1>(1800/60));
        yp1_clist = [yp1_nanone; yp1_over];

        xp1_success = xp1; xp1_success(yp1_clist) = [];
        yp1_success = yp1; yp1_success(yp1_clist) = [];

        xp1_fail = xp1(yp1_clist);
        yp1_fail = yp1(yp1_clist);
        sav_faildat{ilaser} = xp1_fail;

        % regression:
        pev{icond,ilaser} = polyfit(yp1_success,xp1_success,1);
        [rho(icond,ilaser),pval(icond,ilaser)] = corr(yp1_success,xp1_success);

        h{icond,ilaser} = plot(yp1_success,xp1_success,'o'); hold on
        if ilaser == 1
            h{icond,ilaser}.Color = [0 0 0];
            h{icond,ilaser}.MarkerFaceColor = [0 0 0];
            h{icond,ilaser}.MarkerSize = 3;
        elseif ilaser == 2
            h{icond,ilaser}.Color = [0 0.75 0.5];
            h{icond,ilaser}.MarkerFaceColor = [0 0.75 0.5];
            h{icond,ilaser}.MarkerSize = 3;
        end
        
        ax1 = gca;
        xlimcmp = [xlimcmp; ax1.XLim];
        ylimcmp = [ylimcmp; ax1.YLim];

        sub2 = subplot(2,12,icond*3);
        if ilaser == 1
            cla
        end
        if ilaser == 2

            hbox = daboxplot(sav_faildat,'scatter',2,'whiskers',0,'boxalpha',0.7);
            ax2 = gca;
            ax2.XTickLabel = {'Loff','Lon'};

            hbox.bx(1).FaceColor = [0 0 0];
            hbox.bx(1).FaceAlpha = 0.3;

            hbox.bx(2).FaceColor = [0 0.75 0.5];
            hbox.bx(2).FaceAlpha = 0.3;
            hbox.sc(2).MarkerEdgeColor = [1 1 1];
            hbox.sc(2).MarkerFaceColor = [0 0.75 0.5];
        end
        

    end
end
minxlim = min(xlimcmp(:,1));
maxxlim = max(xlimcmp(:,2));
minylim = min(ylimcmp(:,1));
maxylim = max(ylimcmp(:,2));
plotx = minxlim:(maxxlim-minxlim)/100:maxxlim;
for icond = 1:8
    sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
    xlim([minxlim maxxlim])
    ylim([minylim maxylim])

    for ilaser = 1:2
        vala = pev{icond,ilaser}(1);
        valint = pev{icond,ilaser}(2);
        ploty = plotx*vala + valint;
        ptmp = plot(plotx,ploty); hold on
        if ilaser == 1
            ptmp.Color = [0 0 0];
            % ptmp.LineWidth = 2;
            ptmp.LineStyle = ':'
        elseif ilaser == 2
            ptmp.Color = [0 0.75 0.5];
            ptmp.LineStyle = '--';
        end
    end

    title({cond_list{icond},['r=' num2str(rho(icond,1)) ', p=' num2str(pval(icond,1))],['r=' num2str(rho(icond,2)) ', p=' num2str(pval(icond,2))]})


    if icond == 5 || icond == 6 || icond == 7 || icond == 8
        xlabel('capt time (s)')
    end

    if icond == 1 || icond == 5
        ylabel('tot app time (s)')
    end


    sub2 = subplot(2,12,icond*3);
    ylim([minylim maxylim])
end
F = getframe(f);
imwrite(F.cdata, 'scatter_tot_app_time.png', 'png')
close






% % of intercept (time)
f = figure(4437); clf
f.Position = [114 227 1795 625];
h = []; xlimcmp = []; ylimcmp = [];
rho = []; pval = [];
pev = [];
for icond = 1:8
    for ilaser = 1:2
        pev{icond,ilaser} = [];
    end
end
for icond = 1:8
    sav_faildat = [];
    for ilaser = 1:2

        sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
        if ilaser == 1
            cla
        end

        xp1 = int_perc{icond,ilaser};
        yp1 = capt_time{icond,ilaser}/60;
        yp1_nan = isnan(yp1);
        yp1_nanone = find(yp1_nan == 1);
        yp1_over = find(yp1>(1800/60));
        yp1_clist = [yp1_nanone; yp1_over];

        xp1_success = xp1; xp1_success(yp1_clist) = [];
        yp1_success = yp1; yp1_success(yp1_clist) = [];

        xp1_fail = xp1(yp1_clist);
        yp1_fail = yp1(yp1_clist);
        sav_faildat{ilaser} = xp1_fail;

        % regression:
        pev{icond,ilaser} = polyfit(yp1_success,xp1_success,1);
        [rho(icond,ilaser),pval(icond,ilaser)] = corr(yp1_success,xp1_success);

        h{icond,ilaser} = plot(yp1_success,xp1_success,'o'); hold on
        if ilaser == 1
            h{icond,ilaser}.Color = [0 0 0];
            h{icond,ilaser}.MarkerFaceColor = [0 0 0];
            h{icond,ilaser}.MarkerSize = 3;
        elseif ilaser == 2
            h{icond,ilaser}.Color = [0 0.75 0.5];
            h{icond,ilaser}.MarkerFaceColor = [0 0.75 0.5];
            h{icond,ilaser}.MarkerSize = 3;
        end
        
        ax1 = gca;
        xlimcmp = [xlimcmp; ax1.XLim];
        ylimcmp = [ylimcmp; ax1.YLim];

        sub2 = subplot(2,12,icond*3);
        if ilaser == 1
            cla
        end
        if ilaser == 2

            hbox = daboxplot(sav_faildat,'scatter',2,'whiskers',0,'boxalpha',0.7);
            ax2 = gca;
            ax2.XTickLabel = {'Loff','Lon'};

            hbox.bx(1).FaceColor = [0 0 0];
            hbox.bx(1).FaceAlpha = 0.3;

            hbox.bx(2).FaceColor = [0 0.75 0.5];
            hbox.bx(2).FaceAlpha = 0.3;
            hbox.sc(2).MarkerEdgeColor = [1 1 1];
            hbox.sc(2).MarkerFaceColor = [0 0.75 0.5];
        end
        

    end
end
minxlim = min(xlimcmp(:,1));
maxxlim = max(xlimcmp(:,2));
minylim = min(ylimcmp(:,1));
maxylim = max(ylimcmp(:,2));
plotx = minxlim:(maxxlim-minxlim)/100:maxxlim;
for icond = 1:8
    sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
    xlim([minxlim maxxlim])
    ylim([minylim maxylim])

    for ilaser = 1:2
        vala = pev{icond,ilaser}(1);
        valint = pev{icond,ilaser}(2);
        ploty = plotx*vala + valint;
        ptmp = plot(plotx,ploty); hold on
        if ilaser == 1
            ptmp.Color = [0 0 0];
            % ptmp.LineWidth = 2;
            ptmp.LineStyle = ':'
        elseif ilaser == 2
            ptmp.Color = [0 0.75 0.5];
            ptmp.LineStyle = '--';
        end
    end

    title({cond_list{icond},['r=' num2str(rho(icond,1)) ', p=' num2str(pval(icond,1))],['r=' num2str(rho(icond,2)) ', p=' num2str(pval(icond,2))]})


    if icond == 5 || icond == 6 || icond == 7 || icond == 8
        xlabel('capt time (s)')
    end

    if icond == 1 || icond == 5
        ylabel('int % of time')
    end


    sub2 = subplot(2,12,icond*3);
    ylim([minylim maxylim])
end
F = getframe(f);
imwrite(F.cdata, 'scatter_int_perc.png', 'png')
close






% num of int
f = figure(6655); clf
f.Position = [114 227 1795 625];
h = []; xlimcmp = []; ylimcmp = [];
rho = []; pval = [];
pev = [];
for icond = 1:8
    for ilaser = 1:2
        pev{icond,ilaser} = [];
    end
end
for icond = 1:8
    sav_faildat = [];
    for ilaser = 1:2

        sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
        if ilaser == 1
            cla
        end

        xp1 = int_num{icond,ilaser};
        yp1 = capt_time{icond,ilaser}/60;
        yp1_nan = isnan(yp1);
        yp1_nanone = find(yp1_nan == 1);
        yp1_over = find(yp1>(1800/60));
        yp1_clist = [yp1_nanone; yp1_over];

        xp1_success = xp1; xp1_success(yp1_clist) = [];
        yp1_success = yp1; yp1_success(yp1_clist) = [];

        xp1_fail = xp1(yp1_clist);
        yp1_fail = yp1(yp1_clist);
        sav_faildat{ilaser} = xp1_fail;

        % regression:
        pev{icond,ilaser} = polyfit(yp1_success,xp1_success,1);
        [rho(icond,ilaser),pval(icond,ilaser)] = corr(yp1_success,xp1_success);

        h{icond,ilaser} = plot(yp1_success,xp1_success,'o'); hold on
        if ilaser == 1
            h{icond,ilaser}.Color = [0 0 0];
            h{icond,ilaser}.MarkerFaceColor = [0 0 0];
            h{icond,ilaser}.MarkerSize = 3;
        elseif ilaser == 2
            h{icond,ilaser}.Color = [0 0.75 0.5];
            h{icond,ilaser}.MarkerFaceColor = [0 0.75 0.5];
            h{icond,ilaser}.MarkerSize = 3;
        end
        
        ax1 = gca;
        xlimcmp = [xlimcmp; ax1.XLim];
        ylimcmp = [ylimcmp; ax1.YLim];

        sub2 = subplot(2,12,icond*3);
        if ilaser == 1
            cla
        end
        if ilaser == 2

            hbox = daboxplot(sav_faildat,'scatter',2,'whiskers',0,'boxalpha',0.7);
            ax2 = gca;
            ax2.XTickLabel = {'Loff','Lon'};

            hbox.bx(1).FaceColor = [0 0 0];
            hbox.bx(1).FaceAlpha = 0.3;

            hbox.bx(2).FaceColor = [0 0.75 0.5];
            hbox.bx(2).FaceAlpha = 0.3;
            hbox.sc(2).MarkerEdgeColor = [1 1 1];
            hbox.sc(2).MarkerFaceColor = [0 0.75 0.5];
        end
        

    end
end
minxlim = min(xlimcmp(:,1));
maxxlim = max(xlimcmp(:,2));
minylim = min(ylimcmp(:,1));
maxylim = max(ylimcmp(:,2));
plotx = minxlim:(maxxlim-minxlim)/100:maxxlim;
for icond = 1:8
    sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
    xlim([minxlim maxxlim])
    ylim([minylim maxylim])

    for ilaser = 1:2
        vala = pev{icond,ilaser}(1);
        valint = pev{icond,ilaser}(2);
        ploty = plotx*vala + valint;
        ptmp = plot(plotx,ploty); hold on
        if ilaser == 1
            ptmp.Color = [0 0 0];
            % ptmp.LineWidth = 2;
            ptmp.LineStyle = ':'
        elseif ilaser == 2
            ptmp.Color = [0 0.75 0.5];
            ptmp.LineStyle = '--';
        end
    end

    title({cond_list{icond},['r=' num2str(rho(icond,1)) ', p=' num2str(pval(icond,1))],['r=' num2str(rho(icond,2)) ', p=' num2str(pval(icond,2))]})


    if icond == 5 || icond == 6 || icond == 7 || icond == 8
        xlabel('capt time (s)')
    end

    if icond == 1 || icond == 5
        ylabel('int num (#)')
    end


    sub2 = subplot(2,12,icond*3);
    ylim([minylim maxylim])
end
F = getframe(f);
imwrite(F.cdata, 'scatter_int_num.png', 'png')
close



%%

%
% this section compiles all the approach labeled frames and gathers the
% spd, az, and dcrkt - and averages per trial
tmspd = [];
taz = [];
tdcrkt = [];
tcspd = [];
for icond = 1:8
    for ilaser = 1:2
        tmspd{icond,ilaser} = [];
        taz{icond,ilaser} = [];
        tdcrkt{icond,ilaser} = [];
        tcspd{icond,ilaser} = [];
    end
end

% icond = 1;
% targetcond = cond_list{icond};
% targetlaser = 'n';
targetexpctrl = 'Exp';

% first, histograms of lengths of approach distributions
for ifile = 1:size(app_h,2)

    % note that for 'app_d_v8.mat', app_h has the correct index, everything
    % else is +1
    fname = app_h{ifile}.fname;
    fdate = app_h{ifile}.fdate; % only used to locate target_row in df_meta
    fcond = app_h{ifile}.fcond;
    flaser = app_h{ifile}.flaser;
    ftrial = app_h{ifile}.ftrial;

    % if strcmp(fcond,targetcond) == 1 && strcmp(flaser,targetlaser) == 1 

    cond_loc = [];
    for icond = 1:8
        if strcmp(fcond,cond_list{icond})
            cond_loc = icond;
        end
    end
    if isempty(cond_loc)
        disp(['Could not find the correct condition for this trial, something is wrong: ifile #' num2str(ifile)])
    end

    laser_loc = [];
    if strcmp(flaser,'n')
        laser_loc = 1;
    elseif strcmp(flaser,'y')
        laser_loc = 2;
    end

    % get capture time (if any)
    capt_binary = [];
    capt_t = [];
    % find df_meta correspondence
    df_loc = [];
    for idf = 2:size(df_meta,1)
        df_name = df_meta{idf,17};
        df_date = num2str(df_meta{idf,18});
        if length(df_date) == 5
            df_date = ['0' df_date];
        end
        df_cond = df_meta{idf,9};
        df_trial = num2str(df_meta{idf,16});
        df_laser = df_meta{idf,13};
        if strcmp(df_laser,'True')
            df_laser = 'y';
        elseif strcmp(df_laser,'False')
            df_laser = 'n';
        else
            disp('What else is there here other than laser-on vs off???')
        end

        if strcmp(df_name,fname) && strcmp(df_date,fdate) && strcmp(df_cond,fcond) && strcmp(df_trial,ftrial) && strcmp(df_laser,flaser)
            df_loc = idf;
            expctrl = df_meta{df_loc,10};
        end
    end

    if strcmp(expctrl,targetexpctrl)


        disp(['Currently: ' num2str(ifile) ' of ' num2str(size(app_h,2)) ' ... ' fname ', ' fdate ', ' fcond ', ' flaser ', ' ftrial])

        
        file_id = ifile + 0; % correct for 'app_d_v8' here by +1


        for iparam = 1:4 % 4 param readouts (spd, az, dist, cspd)
            tmp_param = app_p_nan{file_id}(iparam,:);
            nan_list = isnan(tmp_param);
            nan_list = find(nan_list == 1);
            tmp_param(nan_list) = [];
            tmp = mean(tmp_param);

            if iparam == 1
                tmspd{cond_loc,laser_loc} = [tmspd{cond_loc,laser_loc}; tmp];
            elseif iparam == 2
                taz{cond_loc,laser_loc} = [taz{cond_loc,laser_loc}; tmp];
            elseif iparam == 3
                tdcrkt{cond_loc,laser_loc} = [tdcrkt{cond_loc,laser_loc}; tmp];
            elseif iparam == 4
                tcspd{cond_loc,laser_loc} = [tcspd{cond_loc,laser_loc}; tmp];
            end

        end


    end
end




% spd
f = figure(1001); clf
f.Position = [114 227 1795 625];
h = []; xlimcmp = []; ylimcmp = [];
rho = []; pval = [];
pev = [];
for icond = 1:8
    for ilaser = 1:2
        pev{icond,ilaser} = [];
    end
end
for icond = 1:8
    sav_faildat = [];
    for ilaser = 1:2

        sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
        if ilaser == 1
            cla
        end

        xp1 = tmspd{icond,ilaser};
        yp1 = capt_time{icond,ilaser}/60;
        yp1_nan = isnan(yp1);
        yp1_nanone = find(yp1_nan == 1);
        yp1_over = find(yp1>(1800/60));
        yp1_clist = [yp1_nanone; yp1_over];

        xp1_success = xp1; xp1_success(yp1_clist) = [];
        yp1_success = yp1; yp1_success(yp1_clist) = [];

        xp1_fail = xp1(yp1_clist);
        yp1_fail = yp1(yp1_clist);
        sav_faildat{ilaser} = xp1_fail;

        % regression:
        pev{icond,ilaser} = polyfit(yp1_success,xp1_success,1);
        [rho(icond,ilaser),pval(icond,ilaser)] = corr(yp1_success,xp1_success);

        h{icond,ilaser} = plot(yp1_success,xp1_success,'o'); hold on
        if ilaser == 1
            h{icond,ilaser}.Color = [0 0 0];
            h{icond,ilaser}.MarkerFaceColor = [0 0 0];
            h{icond,ilaser}.MarkerSize = 3;
        elseif ilaser == 2
            h{icond,ilaser}.Color = [0 0.75 0.5];
            h{icond,ilaser}.MarkerFaceColor = [0 0.75 0.5];
            h{icond,ilaser}.MarkerSize = 3;
        end
        
        ax1 = gca;
        xlimcmp = [xlimcmp; ax1.XLim];
        ylimcmp = [ylimcmp; ax1.YLim];

        sub2 = subplot(2,12,icond*3);
        if ilaser == 1
            cla
        end
        if ilaser == 2

            hbox = daboxplot(sav_faildat,'scatter',2,'whiskers',0,'boxalpha',0.7);
            ax2 = gca;
            ax2.XTickLabel = {'Loff','Lon'};

            hbox.bx(1).FaceColor = [0 0 0];
            hbox.bx(1).FaceAlpha = 0.3;

            hbox.bx(2).FaceColor = [0 0.75 0.5];
            hbox.bx(2).FaceAlpha = 0.3;
            hbox.sc(2).MarkerEdgeColor = [1 1 1];
            hbox.sc(2).MarkerFaceColor = [0 0.75 0.5];
        end
        

    end
end
minxlim = min(xlimcmp(:,1));
maxxlim = max(xlimcmp(:,2));
minylim = min(ylimcmp(:,1));
maxylim = max(ylimcmp(:,2));
plotx = minxlim:(maxxlim-minxlim)/100:maxxlim;
for icond = 1:8
    sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
    xlim([minxlim maxxlim])
    ylim([minylim maxylim])

    for ilaser = 1:2
        vala = pev{icond,ilaser}(1);
        valint = pev{icond,ilaser}(2);
        ploty = plotx*vala + valint;
        ptmp = plot(plotx,ploty); hold on
        if ilaser == 1
            ptmp.Color = [0 0 0];
            % ptmp.LineWidth = 2;
            ptmp.LineStyle = ':'
        elseif ilaser == 2
            ptmp.Color = [0 0.75 0.5];
            ptmp.LineStyle = '--';
        end
    end

    title({cond_list{icond},['r=' num2str(rho(icond,1)) ', p=' num2str(pval(icond,1))],['r=' num2str(rho(icond,2)) ', p=' num2str(pval(icond,2))]})


    if icond == 5 || icond == 6 || icond == 7 || icond == 8
        xlabel('capt time (s)')
    end

    if icond == 1 || icond == 5
        ylabel('avg spd (cm/s)')
    end


    sub2 = subplot(2,12,icond*3);
    ylim([minylim maxylim])
end
F = getframe(f);
imwrite(F.cdata, 'scatter_avg_spd.png', 'png')
close



% angle/az
f = figure(1002); clf
f.Position = [114 227 1795 625];
h = []; xlimcmp = []; ylimcmp = [];
rho = []; pval = [];
pev = [];
for icond = 1:8
    for ilaser = 1:2
        pev{icond,ilaser} = [];
    end
end
for icond = 1:8
    sav_faildat = [];
    for ilaser = 1:2

        sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
        if ilaser == 1
            cla
        end

        xp1 = taz{icond,ilaser};
        yp1 = capt_time{icond,ilaser}/60;
        yp1_nan = isnan(yp1);
        yp1_nanone = find(yp1_nan == 1);
        yp1_over = find(yp1>(1800/60));
        yp1_clist = [yp1_nanone; yp1_over];

        xp1_success = xp1; xp1_success(yp1_clist) = [];
        yp1_success = yp1; yp1_success(yp1_clist) = [];

        xp1_fail = xp1(yp1_clist);
        yp1_fail = yp1(yp1_clist);
        sav_faildat{ilaser} = xp1_fail;

        % regression:
        pev{icond,ilaser} = polyfit(yp1_success,xp1_success,1);
        [rho(icond,ilaser),pval(icond,ilaser)] = corr(yp1_success,xp1_success);

        h{icond,ilaser} = plot(yp1_success,xp1_success,'o'); hold on
        if ilaser == 1
            h{icond,ilaser}.Color = [0 0 0];
            h{icond,ilaser}.MarkerFaceColor = [0 0 0];
            h{icond,ilaser}.MarkerSize = 3;
        elseif ilaser == 2
            h{icond,ilaser}.Color = [0 0.75 0.5];
            h{icond,ilaser}.MarkerFaceColor = [0 0.75 0.5];
            h{icond,ilaser}.MarkerSize = 3;
        end
        
        ax1 = gca;
        xlimcmp = [xlimcmp; ax1.XLim];
        ylimcmp = [ylimcmp; ax1.YLim];

        sub2 = subplot(2,12,icond*3);
        if ilaser == 1
            cla
        end
        if ilaser == 2

            hbox = daboxplot(sav_faildat,'scatter',2,'whiskers',0,'boxalpha',0.7);
            ax2 = gca;
            ax2.XTickLabel = {'Loff','Lon'};

            hbox.bx(1).FaceColor = [0 0 0];
            hbox.bx(1).FaceAlpha = 0.3;

            hbox.bx(2).FaceColor = [0 0.75 0.5];
            hbox.bx(2).FaceAlpha = 0.3;
            hbox.sc(2).MarkerEdgeColor = [1 1 1];
            hbox.sc(2).MarkerFaceColor = [0 0.75 0.5];
        end
        

    end
end
minxlim = min(xlimcmp(:,1));
maxxlim = max(xlimcmp(:,2));
minylim = min(ylimcmp(:,1));
maxylim = max(ylimcmp(:,2));
plotx = minxlim:(maxxlim-minxlim)/100:maxxlim;
for icond = 1:8
    sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
    xlim([minxlim maxxlim])
    ylim([minylim maxylim])

    for ilaser = 1:2
        vala = pev{icond,ilaser}(1);
        valint = pev{icond,ilaser}(2);
        ploty = plotx*vala + valint;
        ptmp = plot(plotx,ploty); hold on
        if ilaser == 1
            ptmp.Color = [0 0 0];
            % ptmp.LineWidth = 2;
            ptmp.LineStyle = ':'
        elseif ilaser == 2
            ptmp.Color = [0 0.75 0.5];
            ptmp.LineStyle = '--';
        end
    end

    title({cond_list{icond},['r=' num2str(rho(icond,1)) ', p=' num2str(pval(icond,1))],['r=' num2str(rho(icond,2)) ', p=' num2str(pval(icond,2))]})


    if icond == 5 || icond == 6 || icond == 7 || icond == 8
        xlabel('capt time (s)')
    end

    if icond == 1 || icond == 5
        ylabel('angle (deg)')
    end


    sub2 = subplot(2,12,icond*3);
    ylim([minylim maxylim])
end
F = getframe(f);
imwrite(F.cdata, 'scatter_avg_az.png', 'png')
close




% dcrkt
f = figure(1003); clf
f.Position = [114 227 1795 625];
h = []; xlimcmp = []; ylimcmp = [];
rho = []; pval = [];
pev = [];
for icond = 1:8
    for ilaser = 1:2
        pev{icond,ilaser} = [];
    end
end
for icond = 1:8
    sav_faildat = [];
    for ilaser = 1:2

        sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
        if ilaser == 1
            cla
        end

        xp1 = tdcrkt{icond,ilaser};
        yp1 = capt_time{icond,ilaser}/60;
        yp1_nan = isnan(yp1);
        yp1_nanone = find(yp1_nan == 1);
        yp1_over = find(yp1>(1800/60));
        yp1_clist = [yp1_nanone; yp1_over];

        xp1_success = xp1; xp1_success(yp1_clist) = [];
        yp1_success = yp1; yp1_success(yp1_clist) = [];

        xp1_fail = xp1(yp1_clist);
        yp1_fail = yp1(yp1_clist);
        sav_faildat{ilaser} = xp1_fail;

        % regression:
        pev{icond,ilaser} = polyfit(yp1_success,xp1_success,1);
        [rho(icond,ilaser),pval(icond,ilaser)] = corr(yp1_success,xp1_success);

        h{icond,ilaser} = plot(yp1_success,xp1_success,'o'); hold on
        if ilaser == 1
            h{icond,ilaser}.Color = [0 0 0];
            h{icond,ilaser}.MarkerFaceColor = [0 0 0];
            h{icond,ilaser}.MarkerSize = 3;
        elseif ilaser == 2
            h{icond,ilaser}.Color = [0 0.75 0.5];
            h{icond,ilaser}.MarkerFaceColor = [0 0.75 0.5];
            h{icond,ilaser}.MarkerSize = 3;
        end
        
        ax1 = gca;
        xlimcmp = [xlimcmp; ax1.XLim];
        ylimcmp = [ylimcmp; ax1.YLim];

        sub2 = subplot(2,12,icond*3);
        if ilaser == 1
            cla
        end
        if ilaser == 2

            hbox = daboxplot(sav_faildat,'scatter',2,'whiskers',0,'boxalpha',0.7);
            ax2 = gca;
            ax2.XTickLabel = {'Loff','Lon'};

            hbox.bx(1).FaceColor = [0 0 0];
            hbox.bx(1).FaceAlpha = 0.3;

            hbox.bx(2).FaceColor = [0 0.75 0.5];
            hbox.bx(2).FaceAlpha = 0.3;
            hbox.sc(2).MarkerEdgeColor = [1 1 1];
            hbox.sc(2).MarkerFaceColor = [0 0.75 0.5];
        end
        

    end
end
minxlim = min(xlimcmp(:,1));
maxxlim = max(xlimcmp(:,2));
minylim = min(ylimcmp(:,1));
maxylim = max(ylimcmp(:,2));
plotx = minxlim:(maxxlim-minxlim)/100:maxxlim;
for icond = 1:8
    sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
    xlim([minxlim maxxlim])
    ylim([minylim maxylim])

    for ilaser = 1:2
        vala = pev{icond,ilaser}(1);
        valint = pev{icond,ilaser}(2);
        ploty = plotx*vala + valint;
        ptmp = plot(plotx,ploty); hold on
        if ilaser == 1
            ptmp.Color = [0 0 0];
            % ptmp.LineWidth = 2;
            ptmp.LineStyle = ':'
        elseif ilaser == 2
            ptmp.Color = [0 0.75 0.5];
            ptmp.LineStyle = '--';
        end
    end

    title({cond_list{icond},['r=' num2str(rho(icond,1)) ', p=' num2str(pval(icond,1))],['r=' num2str(rho(icond,2)) ', p=' num2str(pval(icond,2))]})


    if icond == 5 || icond == 6 || icond == 7 || icond == 8
        xlabel('capt time (s)')
    end

    if icond == 1 || icond == 5
        ylabel('dist (cm)')
    end


    sub2 = subplot(2,12,icond*3);
    ylim([minylim maxylim])
end
F = getframe(f);
imwrite(F.cdata, 'scatter_avg_dcrkt.png', 'png')
close






% crktspd
f = figure(1004); clf
f.Position = [114 227 1795 625];
h = []; xlimcmp = []; ylimcmp = [];
rho = []; pval = [];
pev = [];
for icond = 1:8
    for ilaser = 1:2
        pev{icond,ilaser} = [];
    end
end
for icond = 1:8
    sav_faildat = [];
    for ilaser = 1:2

        sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
        if ilaser == 1
            cla
        end

        xp1 = tcspd{icond,ilaser};
        yp1 = capt_time{icond,ilaser}/60;
        yp1_nan = isnan(yp1);
        yp1_nanone = find(yp1_nan == 1);
        yp1_over = find(yp1>(1800/60));
        yp1_clist = [yp1_nanone; yp1_over];

        xp1_success = xp1; xp1_success(yp1_clist) = [];
        yp1_success = yp1; yp1_success(yp1_clist) = [];

        xp1_fail = xp1(yp1_clist);
        yp1_fail = yp1(yp1_clist);
        sav_faildat{ilaser} = xp1_fail;

        % regression:
        pev{icond,ilaser} = polyfit(yp1_success,xp1_success,1);
        [rho(icond,ilaser),pval(icond,ilaser)] = corr(yp1_success,xp1_success);

        h{icond,ilaser} = plot(yp1_success,xp1_success,'o'); hold on
        if ilaser == 1
            h{icond,ilaser}.Color = [0 0 0];
            h{icond,ilaser}.MarkerFaceColor = [0 0 0];
            h{icond,ilaser}.MarkerSize = 3;
        elseif ilaser == 2
            h{icond,ilaser}.Color = [0 0.75 0.5];
            h{icond,ilaser}.MarkerFaceColor = [0 0.75 0.5];
            h{icond,ilaser}.MarkerSize = 3;
        end
        
        ax1 = gca;
        xlimcmp = [xlimcmp; ax1.XLim];
        ylimcmp = [ylimcmp; ax1.YLim];

        sub2 = subplot(2,12,icond*3);
        if ilaser == 1
            cla
        end
        if ilaser == 2

            hbox = daboxplot(sav_faildat,'scatter',2,'whiskers',0,'boxalpha',0.7);
            ax2 = gca;
            ax2.XTickLabel = {'Loff','Lon'};

            hbox.bx(1).FaceColor = [0 0 0];
            hbox.bx(1).FaceAlpha = 0.3;

            hbox.bx(2).FaceColor = [0 0.75 0.5];
            hbox.bx(2).FaceAlpha = 0.3;
            hbox.sc(2).MarkerEdgeColor = [1 1 1];
            hbox.sc(2).MarkerFaceColor = [0 0.75 0.5];
        end
        

    end
end
minxlim = min(xlimcmp(:,1));
maxxlim = max(xlimcmp(:,2));
minylim = min(ylimcmp(:,1));
maxylim = max(ylimcmp(:,2));
plotx = minxlim:(maxxlim-minxlim)/100:maxxlim;
for icond = 1:8
    sub1 = subplot(2,12,[icond*3-2 icond*3-1]);
    xlim([minxlim maxxlim])
    ylim([minylim maxylim])

    for ilaser = 1:2
        vala = pev{icond,ilaser}(1);
        valint = pev{icond,ilaser}(2);
        ploty = plotx*vala + valint;
        ptmp = plot(plotx,ploty); hold on
        if ilaser == 1
            ptmp.Color = [0 0 0];
            % ptmp.LineWidth = 2;
            ptmp.LineStyle = ':'
        elseif ilaser == 2
            ptmp.Color = [0 0.75 0.5];
            ptmp.LineStyle = '--';
        end
    end

    title({cond_list{icond},['r=' num2str(rho(icond,1)) ', p=' num2str(pval(icond,1))],['r=' num2str(rho(icond,2)) ', p=' num2str(pval(icond,2))]})


    if icond == 5 || icond == 6 || icond == 7 || icond == 8
        xlabel('capt time (s)')
    end

    if icond == 1 || icond == 5
        ylabel('crkt spd (cm/s)')
    end


    sub2 = subplot(2,12,icond*3);
    ylim([minylim maxylim])
end
F = getframe(f);
imwrite(F.cdata, 'scatter_avg_crkt_spd.png', 'png')
close






%% time series plots, using parameters
% param: spd, az, dcrkt, c_spd

% this section compiles all the approach labeled frames and gathers the
% spd, az, and dcrkt for the given frame and puts it into 1 array for each
% condition and laser
tmspd = [];
taz = [];
tdcrkt = [];
% tcspd = [];
for icond = 1:8
    for ilaser = 1:2
        tmspd{icond,ilaser} = [];
        taz{icond,ilaser} = [];
        tdcrkt{icond,ilaser} = [];
        % tcspd{icond,ilaser} = [];
    end
end
tmspd_trial = [];
taz_trial = [];
tdcrkt_trial = [];
for icond = 1:8
    for ilaser = 1:2
        tmspd_trial{icond,ilaser} = [];
        taz_trial{icond,ilaser} = [];
        tdcrkt_trial{icond,ilaser} = [];
        % tcspd{icond,ilaser} = [];
    end
end

capt_time2 = []; % another capt_time (previous one doesn't align for 1 of the conditions?!)
for icond = 1:8
    for ilaser = 1:2
        capt_time2{icond,ilaser} = [];
    end
end
% icond = 1;
% targetcond = cond_list{icond};
% targetlaser = 'n';
targetexpctrl = 'Exp';

% first, histograms of lengths of approach distributions
for ifile = 1:size(app_h,2)

    % note that for 'app_d_v8.mat', app_h has the correct index, everything
    % else is +1
    fname = app_h{ifile}.fname;
    fdate = app_h{ifile}.fdate; % only used to locate target_row in df_meta
    fcond = app_h{ifile}.fcond;
    flaser = app_h{ifile}.flaser;
    ftrial = app_h{ifile}.ftrial;

    % if strcmp(fcond,targetcond) == 1 && strcmp(flaser,targetlaser) == 1 

    cond_loc = [];
    for icond = 1:8
        if strcmp(fcond,cond_list{icond})
            cond_loc = icond;
        end
    end
    if isempty(cond_loc)
        disp(['Could not find the correct condition for this trial, something is wrong: ifile #' num2str(ifile)])
    end

    laser_loc = [];
    if strcmp(flaser,'n')
        laser_loc = 1;
    elseif strcmp(flaser,'y')
        laser_loc = 2;
    end

    % get capture time (if any)
    % capt_binary = [];
    % capt_t = [];
    
    % find df_meta correspondence
    df_loc = [];
    capt_t = 'hello'; % it's a text so it gets replaced; if it doesn't, it'll error out
    for idf = 2:size(df_meta,1)
        df_name = df_meta{idf,17};
        df_date = num2str(df_meta{idf,18});
        if length(df_date) == 5
            df_date = ['0' df_date];
        end
        df_cond = df_meta{idf,9};
        df_trial = num2str(df_meta{idf,16});
        df_laser = df_meta{idf,13};
        if strcmp(df_laser,'True')
            df_laser = 'y';
        elseif strcmp(df_laser,'False')
            df_laser = 'n';
        else
            disp('What else is there here other than laser-on vs off???')
        end

        if strcmp(df_name,fname) && strcmp(df_date,fdate) && strcmp(df_cond,fcond) && strcmp(df_trial,ftrial) && strcmp(df_laser,flaser)
            df_loc = idf;
            expctrl = df_meta{df_loc,10};

            capt_t = round(df_meta{df_loc,2}*60);            

        end
    end

    if strcmp(expctrl,targetexpctrl)

        disp(['Currently: ' num2str(ifile) ' of ' num2str(size(app_h,2)) ' ... ' fname ', ' fdate ', ' fcond ', ' flaser ', ' ftrial])

        
        file_id = ifile + 0; % correct for 'app_d_v8' here by +1


        if isempty(app_d{file_id}) == 0

            tomean_spd = [];
            tomean_az = [];
            tomean_distd = [];

            for iapp = 1:size(app_d{file_id},1)
                tmp_spd = app_stat{file_id}{iapp}.speed;
                tmp_az = app_stat{file_id}{iapp}.az;
                tmp_distd = app_stat{file_id}{iapp}.distd';
                fzero = find(tmp_distd < 0);
                tmp_distd(fzero) = [];
                tmp_spd(fzero) = [];
                tmp_az(fzero) = [];

                tomean_spd = [tomean_spd tmp_spd];
                tomean_az = [tomean_az tmp_az];
                tomean_distd = [tomean_distd tmp_distd];

                tmspd{cond_loc,laser_loc} = [tmspd{cond_loc,laser_loc} tmp_spd];
                taz{cond_loc,laser_loc} = [taz{cond_loc,laser_loc} tmp_az];
                tdcrkt{cond_loc,laser_loc} = [tdcrkt{cond_loc,laser_loc} tmp_distd];
            end

            tmspd_trial{cond_loc,laser_loc} = [tmspd_trial{cond_loc,laser_loc}; mean(tomean_spd)];
            taz_trial{cond_loc,laser_loc} = [taz_trial{cond_loc,laser_loc}; mean(tomean_az)];
            tdcrkt_trial{cond_loc,laser_loc} = [tdcrkt_trial{cond_loc,laser_loc}; mean(tomean_distd)];




            if isempty(capt_t) == 1
                capt_time2{cond_loc,laser_loc} = [capt_time2{cond_loc,laser_loc}; NaN];
            else
                capt_time2{cond_loc,laser_loc} = [capt_time2{cond_loc,laser_loc}; capt_t];
            end
        end

    end
end




% density plots
f = figure(6691); clf
f.Position = [214 178 1596 704];
xlimcmp = [];
ylimcmp = [];
climcmp = [];
for icond = 1:8
    for ilaser = 1:2
        if icond < 5
            subplot(2,11,icond*3-(3-ilaser)); cla
        elseif icond > 4
            subplot(2,11,icond*3-(3-ilaser)-1); cla
        end
        densityScatterChart(tmspd{icond,ilaser},tdcrkt{icond,ilaser});

        ax1 = gca;
        xlimcmp = [xlimcmp; ax1.XLim(1) ax1.XLim(2)];
        ylimcmp = [ylimcmp; ax1.YLim(1) ax1.YLim(2)];
        climcmp = [climcmp; ax1.CLim(1) ax1.CLim(2)];
    end
end
xlimmin = min(xlimcmp(:,1)); xlimmax = max(xlimcmp(:,2));
ylimmin = min(ylimcmp(:,1)); ylimmax = max(ylimcmp(:,2));
climmin = min(climcmp(:,1)); climmax = max(climcmp(:,2));
for icond = 1:8
    for ilaser = 1:2
        if icond < 5
            subplot(2,11,icond*3-(3-ilaser));
        elseif icond > 4
            subplot(2,11,icond*3-(3-ilaser)-1);
        end

        ax1 = gca;
        ax1.XLim = [xlimmin xlimmax];
        ax1.YLim = [ylimmin ylimmax];
        ax1.CLim = [climmin climmax];

        colorbar off
    end
end

% cd('F:\prey_capture_lkhd_analysis\figure\2024_11_12')
cd('F:\prey_capture_lkhd_analysis\figure\2025_01_06c\')
% plot indiv cond
for icond = 1:8
    f = figure(6691+icond); clf
    % f.Position = [214 178 1596 704];
    f.Position = [634 142 821 296];
    xlimcmp = [];
    ylimcmp = [];
    climcmp = [];
    for ilaser = 1:2
        subplot(1,2,ilaser); cla
        densityScatterChart(tmspd{icond,ilaser},tdcrkt{icond,ilaser});
    
        ax1 = gca;
        xlimcmp = [xlimcmp; ax1.XLim(1) ax1.XLim(2)];
        ylimcmp = [ylimcmp; ax1.YLim(1) ax1.YLim(2)];
        climcmp = [climcmp; ax1.CLim(1) ax1.CLim(2)];
    end
    xlimmin = min(xlimcmp(:,1)); xlimmax = max(xlimcmp(:,2));
    ylimmin = min(ylimcmp(:,1)); ylimmax = max(ylimcmp(:,2));
    climmin = min(climcmp(:,1)); climmax = max(climcmp(:,2));
    for ilaser = 1:2
        subplot(1,2,ilaser)

        ax1 = gca;
        % ax1.XLim = [xlimmin xlimmax];
        % ax1.YLim = [ylimmin ylimmax];
        ax1.XLim = [10 25];
        ax1.YLim = [0 15];
        ax1.CLim = [climmin climmax];

        colorbar off

        xlabel('speed (cm/s)')
        if ilaser == 1
            ylabel('distance to crkt (cm)')
            title([cond_list{icond} ', laseroff'])
        elseif ilaser == 2
            title([cond_list{icond} ', laseron'])
        end
    end

    F = getframe(f);
    imwrite(F.cdata, ['heatmap_spd_vs_dist_' cond_list{icond} '.png'], 'png')
    close
end




% speed vs ang
for icond = 1:8
    f = figure(7691+icond); clf
    % f.Position = [214 178 1596 704];
    f.Position = [634 142 821 296];
    xlimcmp = [];
    ylimcmp = [];
    climcmp = [];
    for ilaser = 1:2
        subplot(1,2,ilaser); cla
        densityScatterChart(tmspd{icond,ilaser},taz{icond,ilaser});
    
        ax1 = gca;
        xlimcmp = [xlimcmp; ax1.XLim(1) ax1.XLim(2)];
        ylimcmp = [ylimcmp; ax1.YLim(1) ax1.YLim(2)];
        climcmp = [climcmp; ax1.CLim(1) ax1.CLim(2)];
    end
    xlimmin = min(xlimcmp(:,1)); xlimmax = max(xlimcmp(:,2));
    ylimmin = min(ylimcmp(:,1)); ylimmax = max(ylimcmp(:,2));
    climmin = min(climcmp(:,1)); climmax = max(climcmp(:,2));
    for ilaser = 1:2
        subplot(1,2,ilaser)

        ax1 = gca;
        % ax1.XLim = [xlimmin xlimmax];
        % ax1.YLim = [ylimmin ylimmax];
        ax1.XLim = [5 35];
        ax1.YLim = [0 45];
        ax1.CLim = [climmin climmax];

        % colorbar off

        xlabel('speed (cm/s)')
        if ilaser == 1
            ylabel('azimuth (deg)')
            title([cond_list{icond} ', laseroff'])
        elseif ilaser == 2
            title([cond_list{icond} ', laseron'])
        end
    end

    F = getframe(f);
    imwrite(F.cdata, ['heatmap_spd_vs_az_' cond_list{icond} '.png'], 'png')
    close
end



% dist vs ang
for icond = 1:8
    f = figure(8691+icond); clf
    % f.Position = [214 178 1596 704];
    f.Position = [634 142 821 296];
    xlimcmp = [];
    ylimcmp = [];
    climcmp = [];
    for ilaser = 1:2
        subplot(1,2,ilaser); cla
        densityScatterChart(tdcrkt{icond,ilaser},taz{icond,ilaser});
    
        ax1 = gca;
        xlimcmp = [xlimcmp; ax1.XLim(1) ax1.XLim(2)];
        ylimcmp = [ylimcmp; ax1.YLim(1) ax1.YLim(2)];
        climcmp = [climcmp; ax1.CLim(1) ax1.CLim(2)];
    end
    xlimmin = min(xlimcmp(:,1)); xlimmax = max(xlimcmp(:,2));
    ylimmin = min(ylimcmp(:,1)); ylimmax = max(ylimcmp(:,2));
    climmin = min(climcmp(:,1)); climmax = max(climcmp(:,2));
    for ilaser = 1:2
        subplot(1,2,ilaser)

        ax1 = gca;
        % ax1.XLim = [xlimmin xlimmax];
        % ax1.YLim = [ylimmin ylimmax];
        ax1.XLim = [0 35];
        ax1.YLim = [0 45];
        ax1.CLim = [climmin climmax];

        % colorbar off

        xlabel('distance (cm)')
        if ilaser == 1
            ylabel('azimuth (deg)')
            title([cond_list{icond} ', laseroff'])
        elseif ilaser == 2
            title([cond_list{icond} ', laseron'])
        end
    end

    F = getframe(f);
    imwrite(F.cdata, ['heatmap_dcrkt_vs_az_' cond_list{icond} '.png'], 'png')
    close
end




% speed histogram (does it run faster?)
f = figure(8311); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(tmspd{icond,ilaser}); hold on % divide by 60, because 60 fr/s
        ax1 = gca;
        if ilaser == 2
            xlimcmp = [xlimcmp; ax1.XLim];
        end

    end
end
% minxlim = min(xlimcmp(:,1));
% maxxlim = max(xlimcmp(:,2));
minxlim = 0;
maxxlim = 50;
v = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2
        edges1 = [minxlim:1:maxxlim];
        h{icond,ilaser} = histogram(tmspd{icond,ilaser},edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.06])
        v{icond,ilaser} = vline(mean(tmspd{icond,ilaser},'omitnan')); hold on
        if ilaser == 1
            h{icond,ilaser}.FaceColor = [0.8 0.25 0];
            v{icond,ilaser}.Color = [0.8 0.25 0];
        else
            h{icond,ilaser}.FaceColor = [0 0.25 0.8];
            v{icond,ilaser}.Color = [0 0.25 0.8];
        end
        v{icond,ilaser}.LineStyle = '-';
        v{icond,ilaser}.LineWidth = 1;

    end

    if find(icond == [1 5])
        ylabel('fraction (%/100)')
    end
    if find(icond == [5 6 7 8])
        xlabel('mouse speed (cm/s)')
    end
    legend('laser-n','laser-y')
    title([cond_list{icond} ', Exp'])
end
F = getframe(f);
imwrite(F.cdata, ['hist_spd.png'], 'png')
close

% az histogram (does it look straight-er?)
f = figure(6642); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(taz{icond,ilaser}); hold on % divide by 60, because 60 fr/s
        ax1 = gca;
        if ilaser == 2
            xlimcmp = [xlimcmp; ax1.XLim];
        end

    end
end
% minxlim = min(xlimcmp(:,1));
% maxxlim = max(xlimcmp(:,2));
minxlim = 0;
maxxlim = 45;
v = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2
        edges1 = [minxlim:1:maxxlim];
        h{icond,ilaser} = histogram(taz{icond,ilaser},edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.035])
        v{icond,ilaser} = vline(mean(taz{icond,ilaser},'omitnan')); hold on
        if ilaser == 1
            h{icond,ilaser}.FaceColor = [0.8 0.25 0];
            v{icond,ilaser}.Color = [0.8 0.25 0];
        else
            h{icond,ilaser}.FaceColor = [0 0.25 0.8];
            v{icond,ilaser}.Color = [0 0.25 0.8];
        end
        v{icond,ilaser}.LineStyle = '-';
        v{icond,ilaser}.LineWidth = 1;

    end

    if find(icond == [1 5])
        ylabel('fraction (%/100)')
    end
    if find(icond == [5 6 7 8])
        xlabel('head angle (deg)')
    end
    legend('laser-n','laser-y')
    title([cond_list{icond} ', Exp'])
end
F = getframe(f);
imwrite(F.cdata, ['hist_az.png'], 'png')
close


% dcrkt histogram (does it approach from more-afar? more closer?)
f = figure(5534); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(tdcrkt{icond,ilaser}); hold on % divide by 60, because 60 fr/s
        ax1 = gca;
        if ilaser == 2
            xlimcmp = [xlimcmp; ax1.XLim];
        end

    end
end
% minxlim = min(xlimcmp(:,1));
% maxxlim = max(xlimcmp(:,2));
minxlim = 0;
maxxlim = 55;
v = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2
        edges1 = [minxlim:1:maxxlim];
        h{icond,ilaser} = histogram(tdcrkt{icond,ilaser},edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.15])
        v{icond,ilaser} = vline(mean(tdcrkt{icond,ilaser},'omitnan')); hold on
        if ilaser == 1
            h{icond,ilaser}.FaceColor = [0.8 0.25 0];
            v{icond,ilaser}.Color = [0.8 0.25 0];
        else
            h{icond,ilaser}.FaceColor = [0 0.25 0.8];
            v{icond,ilaser}.Color = [0 0.25 0.8];
        end
        v{icond,ilaser}.LineStyle = '-';
        v{icond,ilaser}.LineWidth = 1;

    end

    if find(icond == [1 5])
        ylabel('fraction (%/100)')
    end
    if find(icond == [5 6 7 8])
        xlabel('dist to crkt (cm)')
    end
    legend('laser-n','laser-y')
    title([cond_list{icond} ', Exp'])
end
F = getframe(f);
imwrite(F.cdata, ['hist_distd.png'], 'png')
close

%
%% time series plots, using parameters (same as above, but exempt <4cm intercept zones)
% param: spd, az, dcrkt, c_spd
tmspd = [];
taz = [];
tdcrkt = [];
% tcspd = [];
for icond = 1:8
    for ilaser = 1:2
        tmspd{icond,ilaser} = [];
        taz{icond,ilaser} = [];
        tdcrkt{icond,ilaser} = [];
        % tcspd{icond,ilaser} = [];
    end
end

% icond = 1;
% targetcond = cond_list{icond};
% targetlaser = 'n';
targetexpctrl = 'Exp';

% first, histograms of lengths of approach distributions
for ifile = 1:size(app_h,2)

    % note that for 'app_d_v8.mat', app_h has the correct index, everything
    % else is +1
    fname = app_h{ifile}.fname;
    fdate = app_h{ifile}.fdate; % only used to locate target_row in df_meta
    fcond = app_h{ifile}.fcond;
    flaser = app_h{ifile}.flaser;
    ftrial = app_h{ifile}.ftrial;

    % if strcmp(fcond,targetcond) == 1 && strcmp(flaser,targetlaser) == 1 

    cond_loc = [];
    for icond = 1:8
        if strcmp(fcond,cond_list{icond})
            cond_loc = icond;
        end
    end
    if isempty(cond_loc)
        disp(['Could not find the correct condition for this trial, something is wrong: ifile #' num2str(ifile)])
    end

    laser_loc = [];
    if strcmp(flaser,'n')
        laser_loc = 1;
    elseif strcmp(flaser,'y')
        laser_loc = 2;
    end

    % get capture time (if any)
    capt_binary = [];
    capt_t = [];
    % find df_meta correspondence
    df_loc = [];
    for idf = 2:size(df_meta,1)
        df_name = df_meta{idf,17};
        df_date = num2str(df_meta{idf,18});
        if length(df_date) == 5
            df_date = ['0' df_date];
        end
        df_cond = df_meta{idf,9};
        df_trial = num2str(df_meta{idf,16});
        df_laser = df_meta{idf,13};
        if strcmp(df_laser,'True')
            df_laser = 'y';
        elseif strcmp(df_laser,'False')
            df_laser = 'n';
        else
            disp('What else is there here other than laser-on vs off???')
        end

        if strcmp(df_name,fname) && strcmp(df_date,fdate) && strcmp(df_cond,fcond) && strcmp(df_trial,ftrial) && strcmp(df_laser,flaser)
            df_loc = idf;
            expctrl = df_meta{df_loc,10};
        end
    end

    if strcmp(expctrl,targetexpctrl)


        disp(['Currently: ' num2str(ifile) ' of ' num2str(size(app_h,2)) ' ... ' fname ', ' fdate ', ' fcond ', ' flaser ', ' ftrial])

        
        file_id = ifile + 0; % correct for 'app_d_v8' here by +1


        if isempty(app_d{file_id}) == 0

            for iapp = 1:size(app_d{file_id},1)
                tmp_spd = app_stat{file_id}{iapp}.speed;
                tmp_az = app_stat{file_id}{iapp}.az;
                tmp_distd = app_stat{file_id}{iapp}.distd';
                ffour = find(tmp_distd < 4);

                tmp_distd(ffour) = [];
                tmp_spd(ffour) = [];
                tmp_az(ffour) = [];

                tmspd{cond_loc,laser_loc} = [tmspd{cond_loc,laser_loc} tmp_spd];
                taz{cond_loc,laser_loc} = [taz{cond_loc,laser_loc} tmp_az];
                tdcrkt{cond_loc,laser_loc} = [tdcrkt{cond_loc,laser_loc} tmp_distd];
            end

        end

    end
end

% speed histogram (does it run faster?)
f = figure(8311); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(tmspd{icond,ilaser}); hold on % divide by 60, because 60 fr/s
        ax1 = gca;
        if ilaser == 2
            xlimcmp = [xlimcmp; ax1.XLim];
        end

    end
end
% minxlim = min(xlimcmp(:,1));
% maxxlim = max(xlimcmp(:,2));
minxlim = 0;
maxxlim = 50;
v = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2
        edges1 = [minxlim:1:maxxlim];
        h{icond,ilaser} = histogram(tmspd{icond,ilaser},edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.06])
        v{icond,ilaser} = vline(mean(tmspd{icond,ilaser},'omitnan')); hold on
        if ilaser == 1
            h{icond,ilaser}.FaceColor = [0.8 0.25 0];
            v{icond,ilaser}.Color = [0.8 0.25 0];
        else
            h{icond,ilaser}.FaceColor = [0 0.25 0.8];
            v{icond,ilaser}.Color = [0 0.25 0.8];
        end
        v{icond,ilaser}.LineStyle = '-';
        v{icond,ilaser}.LineWidth = 1;

    end

    if find(icond == [1 5])
        ylabel('fraction (%/100)')
    end
    if find(icond == [5 6 7 8])
        xlabel('mouse speed (cm/s)')
    end
    legend('laser-n','laser-y')
    title([cond_list{icond} ', Exp'])
end


% az histogram (does it look straight-er?)
f = figure(6642); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(taz{icond,ilaser}); hold on % divide by 60, because 60 fr/s
        ax1 = gca;
        if ilaser == 2
            xlimcmp = [xlimcmp; ax1.XLim];
        end

    end
end
% minxlim = min(xlimcmp(:,1));
% maxxlim = max(xlimcmp(:,2));
minxlim = 0;
maxxlim = 45;
v = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2
        edges1 = [minxlim:1:maxxlim];
        h{icond,ilaser} = histogram(taz{icond,ilaser},edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.035])
        v{icond,ilaser} = vline(mean(taz{icond,ilaser},'omitnan')); hold on
        if ilaser == 1
            h{icond,ilaser}.FaceColor = [0.8 0.25 0];
            v{icond,ilaser}.Color = [0.8 0.25 0];
        else
            h{icond,ilaser}.FaceColor = [0 0.25 0.8];
            v{icond,ilaser}.Color = [0 0.25 0.8];
        end
        v{icond,ilaser}.LineStyle = '-';
        v{icond,ilaser}.LineWidth = 1;

    end

    if find(icond == [1 5])
        ylabel('fraction (%/100)')
    end
    if find(icond == [5 6 7 8])
        xlabel('head angle (deg)')
    end
    legend('laser-n','laser-y')
    title([cond_list{icond} ', Exp'])
end


% dcrkt histogram (does it approach from more-afar? more closer?)
f = figure(5534); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(tdcrkt{icond,ilaser}); hold on % divide by 60, because 60 fr/s
        ax1 = gca;
        if ilaser == 2
            xlimcmp = [xlimcmp; ax1.XLim];
        end

    end
end
% minxlim = min(xlimcmp(:,1));
% maxxlim = max(xlimcmp(:,2));
minxlim = 0;
maxxlim = 55;
v = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2
        edges1 = [minxlim:1:maxxlim];
        h{icond,ilaser} = histogram(tdcrkt{icond,ilaser},edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.15])
        v{icond,ilaser} = vline(mean(tdcrkt{icond,ilaser},'omitnan')); hold on
        if ilaser == 1
            h{icond,ilaser}.FaceColor = [0.8 0.25 0];
            v{icond,ilaser}.Color = [0.8 0.25 0];
        else
            h{icond,ilaser}.FaceColor = [0 0.25 0.8];
            v{icond,ilaser}.Color = [0 0.25 0.8];
        end
        v{icond,ilaser}.LineStyle = '-';
        v{icond,ilaser}.LineWidth = 1;

    end

    if find(icond == [1 5])
        ylabel('fraction (%/100)')
    end
    if find(icond == [5 6 7 8])
        xlabel('dist to crkt (cm)')
    end
    legend('laser-n','laser-y')
    title([cond_list{icond} ', Exp'])
end



%% time series plots, using parameters (same as above, but only the very first 3 frames; where does it all start at?)
% param: spd, az, dcrkt, c_spd
tmspd = [];
taz = [];
tdcrkt = [];
% tcspd = [];
for icond = 1:8
    for ilaser = 1:2
        tmspd{icond,ilaser} = [];
        taz{icond,ilaser} = [];
        tdcrkt{icond,ilaser} = [];
        % tcspd{icond,ilaser} = [];
    end
end

% icond = 1;
% targetcond = cond_list{icond};
% targetlaser = 'n';
targetexpctrl = 'Exp';

% first, histograms of lengths of approach distributions
for ifile = 1:size(app_h,2)

    % note that for 'app_d_v8.mat', app_h has the correct index, everything
    % else is +1
    fname = app_h{ifile}.fname;
    fdate = app_h{ifile}.fdate; % only used to locate target_row in df_meta
    fcond = app_h{ifile}.fcond;
    flaser = app_h{ifile}.flaser;
    ftrial = app_h{ifile}.ftrial;

    % if strcmp(fcond,targetcond) == 1 && strcmp(flaser,targetlaser) == 1 

    cond_loc = [];
    for icond = 1:8
        if strcmp(fcond,cond_list{icond})
            cond_loc = icond;
        end
    end
    if isempty(cond_loc)
        disp(['Could not find the correct condition for this trial, something is wrong: ifile #' num2str(ifile)])
    end

    laser_loc = [];
    if strcmp(flaser,'n')
        laser_loc = 1;
    elseif strcmp(flaser,'y')
        laser_loc = 2;
    end

    % get capture time (if any)
    capt_binary = [];
    capt_t = [];
    % find df_meta correspondence
    df_loc = [];
    for idf = 2:size(df_meta,1)
        df_name = df_meta{idf,17};
        df_date = num2str(df_meta{idf,18});
        if length(df_date) == 5
            df_date = ['0' df_date];
        end
        df_cond = df_meta{idf,9};
        df_trial = num2str(df_meta{idf,16});
        df_laser = df_meta{idf,13};
        if strcmp(df_laser,'True')
            df_laser = 'y';
        elseif strcmp(df_laser,'False')
            df_laser = 'n';
        else
            disp('What else is there here other than laser-on vs off???')
        end

        if strcmp(df_name,fname) && strcmp(df_date,fdate) && strcmp(df_cond,fcond) && strcmp(df_trial,ftrial) && strcmp(df_laser,flaser)
            df_loc = idf;
            expctrl = df_meta{df_loc,10};
        end
    end

    if strcmp(expctrl,targetexpctrl)


        disp(['Currently: ' num2str(ifile) ' of ' num2str(size(app_h,2)) ' ... ' fname ', ' fdate ', ' fcond ', ' flaser ', ' ftrial])

        
        file_id = ifile + 0; % correct for 'app_d_v8' here by +1


        if isempty(app_d{file_id}) == 0

            for iapp = 1:size(app_d{file_id},1)
                tmp_spd = app_stat{file_id}{iapp}.speed;
                tmp_az = app_stat{file_id}{iapp}.az;
                tmp_distd = app_stat{file_id}{iapp}.distd';

                tmp_spd = tmp_spd(1:3);
                tmp_az = tmp_az(1:3);
                tmp_distd = tmp_distd(1:3);

                tmspd{cond_loc,laser_loc} = [tmspd{cond_loc,laser_loc} tmp_spd];
                taz{cond_loc,laser_loc} = [taz{cond_loc,laser_loc} tmp_az];
                tdcrkt{cond_loc,laser_loc} = [tdcrkt{cond_loc,laser_loc} tmp_distd];
            end

        end

    end
end



% speed histogram (does it run faster?)
f = figure(8311); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(tmspd{icond,ilaser}); hold on % divide by 60, because 60 fr/s
        ax1 = gca;
        if ilaser == 2
            xlimcmp = [xlimcmp; ax1.XLim];
        end

    end
end
% minxlim = min(xlimcmp(:,1));
% maxxlim = max(xlimcmp(:,2));
minxlim = 0;
maxxlim = 50;
v = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2
        edges1 = [minxlim:1:maxxlim];
        h{icond,ilaser} = histogram(tmspd{icond,ilaser},edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.1])
        v{icond,ilaser} = vline(mean(tmspd{icond,ilaser},'omitnan')); hold on
        if ilaser == 1
            h{icond,ilaser}.FaceColor = [0.8 0.25 0];
            v{icond,ilaser}.Color = [0.8 0.25 0];
        else
            h{icond,ilaser}.FaceColor = [0 0.25 0.8];
            v{icond,ilaser}.Color = [0 0.25 0.8];
        end
        v{icond,ilaser}.LineStyle = '-';
        v{icond,ilaser}.LineWidth = 1;

    end

    if find(icond == [1 5])
        ylabel('fraction (%/100)')
    end
    if find(icond == [5 6 7 8])
        xlabel('mouse speed (cm/s)')
    end
    legend('laser-n','laser-y')
    title([cond_list{icond} ', Exp'])
end


% az histogram (does it look straight-er?)
f = figure(6642); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(taz{icond,ilaser}); hold on % divide by 60, because 60 fr/s
        ax1 = gca;
        if ilaser == 2
            xlimcmp = [xlimcmp; ax1.XLim];
        end

    end
end
% minxlim = min(xlimcmp(:,1));
% maxxlim = max(xlimcmp(:,2));
minxlim = 0;
maxxlim = 45;
v = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2
        edges1 = [minxlim:1:maxxlim];
        h{icond,ilaser} = histogram(taz{icond,ilaser},edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.1])
        v{icond,ilaser} = vline(mean(taz{icond,ilaser},'omitnan')); hold on
        if ilaser == 1
            h{icond,ilaser}.FaceColor = [0.8 0.25 0];
            v{icond,ilaser}.Color = [0.8 0.25 0];
        else
            h{icond,ilaser}.FaceColor = [0 0.25 0.8];
            v{icond,ilaser}.Color = [0 0.25 0.8];
        end
        v{icond,ilaser}.LineStyle = '-';
        v{icond,ilaser}.LineWidth = 1;

    end

    if find(icond == [1 5])
        ylabel('fraction (%/100)')
    end
    if find(icond == [5 6 7 8])
        xlabel('head angle (deg)')
    end
    L1 = legend('laser-n','laser-y');
    L1.Location = 'northwest';
    title([cond_list{icond} ', Exp'])
end


% dcrkt histogram (does it approach from more-afar? more closer?)
f = figure(5534); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(tdcrkt{icond,ilaser}); hold on % divide by 60, because 60 fr/s
        ax1 = gca;
        if ilaser == 2
            xlimcmp = [xlimcmp; ax1.XLim];
        end

    end
end
% minxlim = min(xlimcmp(:,1));
% maxxlim = max(xlimcmp(:,2));
minxlim = 0;
maxxlim = 55;
v = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2
        edges1 = [minxlim:1:maxxlim];
        h{icond,ilaser} = histogram(tdcrkt{icond,ilaser},edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.15])
        v{icond,ilaser} = vline(mean(tdcrkt{icond,ilaser},'omitnan')); hold on
        if ilaser == 1
            h{icond,ilaser}.FaceColor = [0.8 0.25 0];
            v{icond,ilaser}.Color = [0.8 0.25 0];
        else
            h{icond,ilaser}.FaceColor = [0 0.25 0.8];
            v{icond,ilaser}.Color = [0 0.25 0.8];
        end
        v{icond,ilaser}.LineStyle = '-';
        v{icond,ilaser}.LineWidth = 1;

    end

    if find(icond == [1 5])
        ylabel('fraction (%/100)')
    end
    if find(icond == [5 6 7 8])
        xlabel('dist to crkt (cm)')
    end
    legend('laser-n','laser-y')
    title([cond_list{icond} ', Exp'])
end



%% time stationary/perc stationary
tstationary = [];
pstationary = [];
for icond = 1:8
    for ilaser = 1:2
        tstationary{icond,ilaser} = [];
        pstationary{icond,ilaser} = [];
    end
end

% icond = 1;
% targetcond = cond_list{icond};
% targetlaser = 'n';
targetexpctrl = 'Exp';

% first, histograms of lengths of approach distributions
for ifile = 1:size(app_h,2)

    % note that for 'app_d_v8.mat', app_h has the correct index, everything
    % else is +1
    fname = app_h{ifile}.fname;
    fdate = app_h{ifile}.fdate; % only used to locate target_row in df_meta
    fcond = app_h{ifile}.fcond;
    flaser = app_h{ifile}.flaser;
    ftrial = app_h{ifile}.ftrial;

    % if strcmp(fcond,targetcond) == 1 && strcmp(flaser,targetlaser) == 1 

    cond_loc = [];
    for icond = 1:8
        if strcmp(fcond,cond_list{icond})
            cond_loc = icond;
        end
    end
    if isempty(cond_loc)
        disp(['Could not find the correct condition for this trial, something is wrong: ifile #' num2str(ifile)])
    end

    laser_loc = [];
    if strcmp(flaser,'n')
        laser_loc = 1;
    elseif strcmp(flaser,'y')
        laser_loc = 2;
    end

    % get capture time (if any)
    capt_binary = [];
    capt_t = [];
    % find df_meta correspondence
    df_loc = [];
    for idf = 2:size(df_meta,1)
        df_name = df_meta{idf,17};
        df_date = num2str(df_meta{idf,18});
        if length(df_date) == 5
            df_date = ['0' df_date];
        end
        df_cond = df_meta{idf,9};
        df_trial = num2str(df_meta{idf,16});
        df_laser = df_meta{idf,13};
        if strcmp(df_laser,'True')
            df_laser = 'y';
        elseif strcmp(df_laser,'False')
            df_laser = 'n';
        else
            disp('What else is there here other than laser-on vs off???')
        end

        if strcmp(df_name,fname) && strcmp(df_date,fdate) && strcmp(df_cond,fcond) && strcmp(df_trial,ftrial) && strcmp(df_laser,flaser)
            df_loc = idf;
            expctrl = df_meta{df_loc,10};
        end
    end

    if strcmp(expctrl,targetexpctrl)


        disp(['Currently: ' num2str(ifile) ' of ' num2str(size(app_h,2)) ' ... ' fname ', ' fdate ', ' fcond ', ' flaser ', ' ftrial])

        
        file_id = ifile + 0; % correct for 'app_d_v8' here by +1


        if isempty(app_d{file_id}) == 0

            for iapp = 1:size(app_d{file_id},1)
                % tmp_spd = app_stat{file_id}{iapp}.speed;
                tmp_spd = app_p_nan{file_id}(1,:);
                t_stationary = length(find(tmp_spd<4)); % below 4cm/s
                tot_len = length(tmp_spd)-sum(isnan(tmp_spd));
                p_stationary = (t_stationary)/tot_len;

                tstationary{cond_loc,laser_loc} = [tstationary{cond_loc,laser_loc} t_stationary];
                pstationary{cond_loc,laser_loc} = [pstationary{cond_loc,laser_loc} p_stationary];
            end

        end

    end
end

% dcrkt histogram (does it approach from more-afar? more closer?)
f = figure(9921); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(tstationary{icond,ilaser}); hold on % divide by 60, because 60 fr/s
        ax1 = gca;
        if ilaser == 2
            xlimcmp = [xlimcmp; ax1.XLim];
        end

    end
end
minxlim = min(xlimcmp(:,1));
maxxlim = max(xlimcmp(:,2));
v = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2
        edges1 = [minxlim:75:maxxlim];
        h{icond,ilaser} = histogram(tstationary{icond,ilaser},edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.5])
        v{icond,ilaser} = vline(mean(tstationary{icond,ilaser},'omitnan')); hold on
        if ilaser == 1
            h{icond,ilaser}.FaceColor = [0.8 0.25 0];
            v{icond,ilaser}.Color = [0.8 0.25 0];
        else
            h{icond,ilaser}.FaceColor = [0 0.25 0.8];
            v{icond,ilaser}.Color = [0 0.25 0.8];
        end
        v{icond,ilaser}.LineStyle = '-';
        v{icond,ilaser}.LineWidth = 1;

    end

    if find(icond == [1 5])
        ylabel('fraction (%/100)')
    end
    if find(icond == [5 6 7 8])
        xlabel('time stationary (fr)')
    end
    legend('laser-n','laser-y')
    title([cond_list{icond} ', Exp'])
end


% 
f = figure(9221); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(pstationary{icond,ilaser}); hold on % divide by 60, because 60 fr/s
        ax1 = gca;
        if ilaser == 2
            xlimcmp = [xlimcmp; ax1.XLim];
        end

    end
end
minxlim = min(xlimcmp(:,1));
maxxlim = max(xlimcmp(:,2));
v = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2
        edges1 = [minxlim:0.025:maxxlim];
        h{icond,ilaser} = histogram(pstationary{icond,ilaser},edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.5])
        v{icond,ilaser} = vline(mean(pstationary{icond,ilaser},'omitnan')); hold on
        if ilaser == 1
            h{icond,ilaser}.FaceColor = [0.8 0.25 0];
            v{icond,ilaser}.Color = [0.8 0.25 0];
        else
            h{icond,ilaser}.FaceColor = [0 0.25 0.8];
            v{icond,ilaser}.Color = [0 0.25 0.8];
        end
        v{icond,ilaser}.LineStyle = '-';
        v{icond,ilaser}.LineWidth = 1;

    end

    if find(icond == [1 5])
        ylabel('fraction (%/100)')
    end
    if find(icond == [5 6 7 8])
        xlabel('time stationary (fr)')
    end
    legend('laser-n','laser-y')
    title([cond_list{icond} ', Exp'])
end

