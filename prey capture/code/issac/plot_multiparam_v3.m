% approach freq, intercept given approach, capture given intercept

addpath('f:\prey_capture_lkhd_analysis\code')  
addpath('f:\code')  
addpath('f:\code\hline_vline\')

% load 'app_d'
% load('f:\prey_capture_lkhd_analysis\d_save\app_d_v9.mat')
% load('f:\prey_capture_lkhd_analysis\d_save\app_d_v11.mat')
load('f:\prey_capture_lkhd_analysis\d_save\app_d_v15.mat')
og_app_int = app_int;
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
int_freq = [];
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
        int_freq{icond,ilaser} = [];
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
% targetexpctrl = 'Exp';
targetexpctrl = 'Ctrl';

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
            vid_len = size(app_p{file_id},2);
            for iapp = 1:size(tmp_appd,1)
                st = tmp_appd(iapp,1);
                ed = tmp_appd(iapp,2);
                ed1 = ed;
                if ed == vid_len
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
                    int_freq{icond,ilaser} = [int_freq{icond,ilaser}; freq_int];
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
        
                    int_perc{icond,ilaser} = [int_perc{icond,ilaser}; 0];
                    int_num{icond,ilaser} = [int_num{icond,ilaser}; 0];
                    int_freq{icond,ilaser} = [int_freq{icond,ilaser}; 0];
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
                end
            end
        end

    end
    % end
end


    
%% Figures

% frequency of approach (count as absolute # divided by time)
f = figure(1829); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(app_freq{icond,ilaser}/60); hold on % divide by 60, because 60 fr/s
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
        edges1 = [minxlim:0.00001:maxxlim];
        h{icond,ilaser} = histogram(app_freq{icond,ilaser}/60,edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.25])
        v{icond,ilaser} = vline(mean(app_freq{icond,ilaser}/60)); hold on
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
        xlabel('approach freq (#/s)')
    end
    legend('laser-n','laser-y')
    title([cond_list{icond} ', Exp'])
end

% percentage of time spent in approach
f = figure(1849); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(app_perc{icond,ilaser}); hold on 
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
        edges1 = [minxlim:0.01:maxxlim];
        h{icond,ilaser} = histogram(app_perc{icond,ilaser},edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.2])
        v{icond,ilaser} = vline(mean(app_perc{icond,ilaser})); hold on
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
        xlabel('approach time (%/100)')
    end
    legend('laser-n','laser-y')
    title([cond_list{icond} ', Exp'])
end

% probability of intercept given approach
f = figure(1949); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(int_per_app{icond,ilaser}); hold on 
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
        edges1 = [minxlim:0.05:maxxlim];
        h{icond,ilaser} = histogram(int_per_app{icond,ilaser},edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.4])
        v{icond,ilaser} = vline(mean(int_per_app{icond,ilaser},'omitnan')); hold on
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
        xlabel('int per approach (%/100)')
    end
    l1 = legend('laser-n','laser-y');
    l1.Location = 'northwest';
    title([cond_list{icond} ', Exp'])
end

% histogram of num of approach
f = figure(1711); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(app_num{icond,ilaser}); hold on % divide by 60, because 60 fr/s
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
        edges1 = [minxlim:1:maxxlim];
        h{icond,ilaser} = histogram(app_num{icond,ilaser},edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.25])
        v{icond,ilaser} = vline(mean(app_num{icond,ilaser})); hold on
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
        xlabel('approach num (#)')
    end
    legend('laser-n','laser-y')
    title([cond_list{icond} ', Exp'])
end

% bar graph of probability of intercept given approach (but calculated for
% entire condition)
% probability of intercept given approach
f = figure(1299); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla

    ilaser = 1;
    tmp1 = sum(int_per_app_stat{icond,ilaser});
    sum11 = tmp1(1);
    sum12 = tmp1(2);

    ilaser = 2;
    tmp2 = sum(int_per_app_stat{icond,ilaser});
    sum21 = tmp2(1);
    sum22 = tmp2(2);

    h{icond} = bar([sum12/sum11 sum22/sum21]);

    ax1 = gca;
    if ilaser == 2
        xlimcmp = [xlimcmp; ax1.XLim];
    end

    ylim([0 0.95])
end
for icond = 1:8
    subplot(2,4,icond)

    ax1 = gca;

    ax1.XTickLabel = {'laser-n','laser-y'};

    if find(icond == [1 5])
        ylabel('fraction (%/100)')
    end
    if find(icond == [5 6 7 8])
        xlabel('int per app (%/100)')
    end
    
    title([cond_list{icond} ', Exp'])
    
end

% probability of intercept per approach (duplicate)
f = figure(1021); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(int_per_app{icond,ilaser}); hold on 
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
        edges1 = [minxlim:0.05:maxxlim];
        h{icond,ilaser} = histogram(int_per_app{icond,ilaser},edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.4])
        v{icond,ilaser} = vline(mean(int_per_app{icond,ilaser},'omitnan')); hold on
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
        xlabel('int per approach (%/100)')
    end
    l1 = legend('laser-n','laser-y');
    l1.Location = 'northwest';
    title([cond_list{icond} ', Exp'])
end

% histogram for capture per intercept (this does not necessarily mean that
% an intercept led to the capture, however!)

% [[[ this is impossible to plot as there are trials where mouse captures
% the cricket without an intercept, meaning the conditions for 'approach'
% wasn't satisfied but the mouse caught the cricket anyway ]]]
% above has been fixed: capt_binary is now only considering prior to 30s
% mark as success

capt_per_int_t = []; % per trial

for icond = 1:8
    for ilaser = 1:2
        capt_per_int_t{icond,ilaser} = [];
    end
end
for icond = 1:8
    for ilaser = 1:2
        for itrial = 1:size(capt_bin{icond,ilaser},1)
            
            tmp_n_capt = capt_bin{icond,ilaser}(itrial);
            tmp_n_int = int_num{icond,ilaser}(itrial);

            if capt_bin{icond,ilaser}(itrial) == 0
                capt_per_int_t{icond,ilaser} = [capt_per_int_t{icond,ilaser}; 0];
            elseif capt_bin{icond,ilaser}(itrial) == 1
                capt_per_int_t{icond,ilaser} = [capt_per_int_t{icond,ilaser}; tmp_n_capt/tmp_n_int];
            else
                disp('What else is there??')
            end
        end
    end
end

% there is Inf in data, meaning capture happened without intercept (or
% more likely, an approach altogether)
f = figure(9283); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(capt_per_int_t{icond,ilaser}); hold on 
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
        edges1 = [minxlim:0.02:maxxlim];
        h{icond,ilaser} = histogram(capt_per_int_t{icond,ilaser},edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.9])
        v{icond,ilaser} = vline(mean(capt_per_int_t{icond,ilaser},'omitnan')); hold on
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
        xlabel('capt per int (%/100)')
    end
    l1 = legend('laser-n','laser-y');
    l1.Location = 'northeast';
    title([cond_list{icond} ', Exp'])
end


% bar graph for the same info, but now, you lose the trial information,
% grouped every trial into condition
capt_per_int = [];
for icond = 1:8
    for ilaser = 1:2
        tmp_n_capt = sum(capt_bin{icond,ilaser});
        tmp_n_int = sum(int_num{icond,ilaser});
        capt_per_int(icond,ilaser) = tmp_n_capt/tmp_n_int;
    end
end


f = figure(1188); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla

    h{icond} = bar(capt_per_int(icond,:));

    ax1 = gca;
    xlimcmp = [xlimcmp; ax1.XLim];

    ylim([0 0.4])
end
for icond = 1:8
    subplot(2,4,icond)

    ax1 = gca;

    ax1.XTickLabel = {'laser-n','laser-y'};

    if find(icond == [1 5])
        ylabel('capture per int (%/100)')
    end
    
    title([cond_list{icond} ', Exp'])
    
end


% lengths of each approach (count as absolute #/time)
f = figure(1219); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(app_len{icond,ilaser}/60); hold on % divide by 60, because 60 fr/s
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
        edges1 = [minxlim:0.1:maxxlim];
        h{icond,ilaser} = histogram(app_len{icond,ilaser}/60,edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.35])
        v{icond,ilaser} = vline(mean(app_len{icond,ilaser}/60,'omitnan')); hold on
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
        xlabel('approach length (s)')
    end
    legend('laser-n','laser-y')
    title([cond_list{icond} ', Exp'])
end



% total time of approach (#/time)
f = figure(3339); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(app_tot{icond,ilaser}/60); hold on % divide by 60, because 60 fr/s
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
        edges1 = [minxlim:0.5:maxxlim];
        h{icond,ilaser} = histogram(app_tot{icond,ilaser}/60,edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.25])
        v{icond,ilaser} = vline(mean(app_tot{icond,ilaser}/60)); hold on
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
        xlabel('approach total time (s)')
    end
    legend('laser-n','laser-y')
    title([cond_list{icond} ', Exp'])
end




% time spent doing intercept (% of total time of trial)
f = figure(2812); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(int_perc{icond,ilaser}); hold on % divide by 60, because 60 fr/s
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
        edges1 = [minxlim:0.01:maxxlim];
        h{icond,ilaser} = histogram(int_perc{icond,ilaser},edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.25])
        v{icond,ilaser} = vline(mean(int_perc{icond,ilaser})); hold on
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
        xlabel('time spent on intercept')
    end
    legend('laser-n','laser-y')
    title([cond_list{icond} ', Exp'])
end


% bar graph of capt binary
f = figure(9938); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
capt_rate_plt = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla

    ilaser = 1;
    tmp1 = sum(capt_bin{icond,ilaser});
    len1 = length(capt_bin{icond,ilaser});

    ilaser = 2;
    tmp2 = sum(capt_bin{icond,ilaser});
    len2 = length(capt_bin{icond,ilaser});

    h{icond} = bar([tmp1/len1 tmp2/len2]);

    capt_rate_plt = [capt_rate_plt; tmp1/len1 tmp2/len2];

    ax1 = gca;
    if ilaser == 2
        xlimcmp = [xlimcmp; ax1.XLim];
    end

    ylim([0 0.95])
end
for icond = 1:8
    subplot(2,4,icond)

    ax1 = gca;

    ax1.XTickLabel = {'laser-n','laser-y'};

    if find(icond == [1 5])
        ylabel('fraction (%/100)')
    end
    if find(icond == [5 6 7 8])
        xlabel('capture rate (%)')
    end
    
    title([cond_list{icond} ', Exp'])
    
end



% time spent doing intercept (% of total time of trial) 
f = figure(2812); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(int_perc{icond,ilaser}); hold on % divide by 60, because 60 fr/s
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
        edges1 = [minxlim:0.01:maxxlim];
        h{icond,ilaser} = histogram(int_perc{icond,ilaser},edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 0.25])
        v{icond,ilaser} = vline(mean(int_perc{icond,ilaser})); hold on
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
        xlabel('time spent on intercept')
    end
    legend('laser-n','laser-y')
    title([cond_list{icond} ', Exp'])
end


% success rate of approach, per app
f = figure(1839); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(success_app{icond,ilaser}); hold on % divide by 60, because 60 fr/s
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
        edges1 = [minxlim:0.05:maxxlim];
        h{icond,ilaser} = histogram(success_app{icond,ilaser},edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 1])
        v{icond,ilaser} = vline(mean(success_app{icond,ilaser},'omitnan')); hold on
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
        xlabel('approach success rate')
    end
    legend('laser-n','laser-y')
    title([cond_list{icond} ', Exp'])
end



% success rate of intercept, per intercept
f = figure(1837); clf
f.Position = [214 178 1596 704];
h = []; xlimcmp = [];
for icond = 1:8
    sub1 = subplot(2,4,icond); cla
    for ilaser = 1:2

        h{icond,ilaser} = histogram(success_int{icond,ilaser}); hold on % divide by 60, because 60 fr/s
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
        edges1 = [minxlim:0.05:maxxlim];
        h{icond,ilaser} = histogram(success_int{icond,ilaser},edges1); hold on
        h{icond,ilaser}.Normalization = 'probability';

        ylim([0 1])
        v{icond,ilaser} = vline(mean(success_int{icond,ilaser},'omitnan')); hold on
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
        xlabel('intercept success rate')
    end
    legend('laser-n','laser-y')
    title([cond_list{icond} ', Exp'])
end

%% time series plots, using parameters
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
                fzero = find(tmp_distd < 0);
                tmp_distd(fzero) = [];

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



%% time series plots, using parameters (same as above, but only the very LAST 3 frames; where does it all start at?)
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

                tmp_spd = tmp_spd(end-2:end);
                tmp_az = tmp_az(end-2:end);
                tmp_distd = tmp_distd(end-2:end);

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
                t_stationary = length(find(tmp_spd<5)); % below 5cm/s
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

