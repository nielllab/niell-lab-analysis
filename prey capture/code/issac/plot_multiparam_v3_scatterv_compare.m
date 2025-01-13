% approach freq, intercept given approach, capture given intercept

addpath('f:\prey_capture_lkhd_analysis\code')  
addpath('f:\code')  
addpath('f:\code\hline_vline\')
addpath('f:\code\povilaskarvelis\daboxplot\')
addpath('F:\code\povilaskarvelis\dabarplot\')
addpath('f:\code\densityScatterChart-1.2.0.0\')

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
int_num = []; % number of intercepts
int_freq = []; % frequency of intercept, #/time
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
        int_freq{icond,ilaser} = [];
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
                    app_freq{icond,ilaser} = [app_freq{icond,ilaser}; 0];
        
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
% bar plots - doesn't discriminate success vs fail trials
% want to answer this:

% - % - % - %

% Mice make less intercept/approach attempts when laser is on, and in
% difficult conditions. True or False?

% - % - % - %

sav_on = 1;
% cd('F:\prey_capture_lkhd_analysis\figure\2024_12_02\') % save here (save folder sav fldr)
cd('F:\prey_capture_lkhd_analysis\figure\2025_01_06c\')
condition_names = cond_list;
group_names = {'laseroff','laseron'};

%%
% approaches:
% 
% freq (app_freq)
% rearrange:
tdat = [];

hd = app_freq;
numd = [];
for icond = 1:8
    for ilaser = 1:2
        numd(icond,ilaser) = length(hd{icond,ilaser});
    end
end
maxn = max(numd(:));
nant = nan(maxn,1);
group_inx = [ones(maxn,1); ones(maxn,1)*2]; % need to label which ones are laser off and laser on
laseroffc = []; % going to join with 'laseronc' to complete the 'tdat' which is the 'data1'
laseronc = [];
for icond = 1:8
    for ilaser = 1:2
        hdt = hd{icond,ilaser};

        hdt_len = length(hdt);

        tmpplug = nant;
        tmpplug(1:hdt_len) = hdt;

        if ilaser == 1
            laseroffc(:,icond) = tmpplug;
        elseif ilaser == 2
            laseronc(:,icond) = tmpplug;
        end
    end
end
tdat = [laseroffc; laseronc];
data1 = tdat;

jj = figure(101); clf
jj.Position = [1026 259 755 459];
h = dabarplot(data1,'groups',group_inx,...
    'xtlabels', condition_names,'errorbars',0,...
    'scatter',1,'scattersize',15,'scatteralpha',0.5,...
    'barspacing',0.8,'legend',group_names); 
ylabel('app freq (#/time)');
% yl = ylim; ylim([yl(1), yl(2)+2]);  % make more space for the legend
set(gca,'FontSize',11);
xlabel('conditions')
title('approach freq (#/times (s))')

if sav_on == 1
    F = getframe(jj);
    imwrite(F.cdata, 'barscatter_app_freq.png', 'png')
    close
end

%%
% approaches:
% 
% perc (app_perc)
% rearrange:
tdat = [];

hd = app_perc;
numd = [];
for icond = 1:8
    for ilaser = 1:2
        numd(icond,ilaser) = length(hd{icond,ilaser});
    end
end
maxn = max(numd(:));
nant = nan(maxn,1);
group_inx = [ones(maxn,1); ones(maxn,1)*2]; % need to label which ones are laser off and laser on
laseroffc = []; % going to join with 'laseronc' to complete the 'tdat' which is the 'data1'
laseronc = [];
for icond = 1:8
    for ilaser = 1:2
        hdt = hd{icond,ilaser};

        hdt_len = length(hdt);

        tmpplug = nant;
        tmpplug(1:hdt_len) = hdt;

        if ilaser == 1
            laseroffc(:,icond) = tmpplug;
        elseif ilaser == 2
            laseronc(:,icond) = tmpplug;
        end
    end
end
tdat = [laseroffc; laseronc];
data1 = tdat;

jj = figure(102); clf
jj.Position = [1026 259 755 459];
h = dabarplot(data1,'groups',group_inx,...
    'xtlabels', condition_names,'errorbars',0,...
    'scatter',1,'scattersize',15,'scatteralpha',0.5,...
    'barspacing',0.8,'legend',group_names); 
ylabel('app perc (%*100)');
% yl = ylim; ylim([yl(1), yl(2)+2]);  % make more space for the legend
set(gca,'FontSize',11);
xlabel('conditions')
title('approach percentage of time (%*100)')

if sav_on == 1
    F = getframe(jj);
    imwrite(F.cdata, 'barscatter_app_perc.png', 'png')
    close
end

%%
% intercepts:
% 
% freq (int_freq)
% rearrange:
tdat = [];

hd = int_freq;
numd = [];
for icond = 1:8
    for ilaser = 1:2
        numd(icond,ilaser) = length(hd{icond,ilaser});
    end
end
maxn = max(numd(:));
nant = nan(maxn,1);
group_inx = [ones(maxn,1); ones(maxn,1)*2]; % need to label which ones are laser off and laser on
laseroffc = []; % going to join with 'laseronc' to complete the 'tdat' which is the 'data1'
laseronc = [];
for icond = 1:8
    for ilaser = 1:2
        hdt = hd{icond,ilaser};

        hdt_len = length(hdt);

        tmpplug = nant;
        tmpplug(1:hdt_len) = hdt;

        if ilaser == 1
            laseroffc(:,icond) = tmpplug;
        elseif ilaser == 2
            laseronc(:,icond) = tmpplug;
        end
    end
end
tdat = [laseroffc; laseronc];
data1 = tdat;

jj = figure(201); clf
jj.Position = [1026 259 755 459];
h = dabarplot(data1,'groups',group_inx,...
    'xtlabels', condition_names,'errorbars',0,...
    'scatter',1,'scattersize',15,'scatteralpha',0.5,...
    'barspacing',0.8,'legend',group_names); 
ylabel('int freq (#/time)');
% yl = ylim; ylim([yl(1), yl(2)+2]);  % make more space for the legend
set(gca,'FontSize',11);
xlabel('conditions')
title('intercept freq (#/times (s))')

if sav_on == 1
    F = getframe(jj);
    imwrite(F.cdata, 'barscatter_int_freq.png', 'png')
    close
end

%%
% intercepts:
% 
% perc (int_perc)
% rearrange:
tdat = [];

hd = int_perc;
numd = [];
for icond = 1:8
    for ilaser = 1:2
        numd(icond,ilaser) = length(hd{icond,ilaser});
    end
end
maxn = max(numd(:));
nant = nan(maxn,1);
group_inx = [ones(maxn,1); ones(maxn,1)*2]; % need to label which ones are laser off and laser on
laseroffc = []; % going to join with 'laseronc' to complete the 'tdat' which is the 'data1'
laseronc = [];
for icond = 1:8
    for ilaser = 1:2
        hdt = hd{icond,ilaser};

        hdt_len = length(hdt);

        tmpplug = nant;
        tmpplug(1:hdt_len) = hdt;

        if ilaser == 1
            laseroffc(:,icond) = tmpplug;
        elseif ilaser == 2
            laseronc(:,icond) = tmpplug;
        end
    end
end
tdat = [laseroffc; laseronc];
data1 = tdat;

jj = figure(202); clf
jj.Position = [1026 259 755 459];
h = dabarplot(data1,'groups',group_inx,...
    'xtlabels', condition_names,'errorbars',0,...
    'scatter',1,'scattersize',15,'scatteralpha',0.5,...
    'barspacing',0.8,'legend',group_names); 
ylabel('int perc (%*100)');
% yl = ylim; ylim([yl(1), yl(2)+2]);  % make more space for the legend
set(gca,'FontSize',11);
xlabel('conditions')
title('intercept percentage of time (%*100)')

if sav_on == 1
    F = getframe(jj);
    imwrite(F.cdata, 'barscatter_int_perc.png', 'png')
    close
end

%%

% i want to know:
% Once the mouse reaches intercept zone, success rate is the same. True or
% false?

% strategy:
% grab multiple trials at once, count the capture number and divide it by
% number of approaches or number of intercepts (hence, this analysis
% ignores the time dimension; the fact that longer trials have more
% approaches/higher time percentage of approaches, etc.)

% Bootstrap multiple trials
bstrap_num = 5; % how many trials
n_pts = 100; % repeat this many

% use 'capt_bin' as it already calculated capture or not (under 1800 fr; laser period)

%%
% for approach:
n_dat = [];
for icond = 1:8
    for ilaser = 1:2
        n_dat{icond,ilaser} = [];
    end
end

for icond = 1:8
    for ilaser = 1:2
        n_trials = length(capt_bin{icond,ilaser});
        app_num_tmp = app_num{icond,ilaser};
        capt_bin_tmp = capt_bin{icond,ilaser};

        for iter = 1:n_pts

            try_i = randi(n_trials,1,bstrap_num);

            tot_app_num = sum(app_num_tmp(try_i));
            tot_capt_bin = sum(capt_bin_tmp(try_i));

            % success/num-app
            s_rate = tot_capt_bin/tot_app_num;

            n_dat{icond,ilaser} = [n_dat{icond,ilaser}; s_rate];

            if rem(iter,100) == 0
                if ilaser == 1
                    laser_txt = 'laser off';
                elseif ilaser == 2
                    laser_txt = 'laser on';
                end
                disp(['cond: ' num2str(icond) ', ilaser: ' laser_txt ', num iter: ' num2str(iter) ' of ' num2str(n_pts)])
            end
        end

    end
end
app_dat = n_dat;

f = figure(3201); clf
f.Position = [854 94 984 830]; % 1 monitor setup
cnt = 0;
h_dat = [];
ymax = []; % used to set the ylim afterwards
med_f = [];
avg_f = [];
for icond = 1:8
    for ilaser = 1:2
        cnt = cnt + 1;
        subplot(8,2,cnt); cla
        h_dat{icond,ilaser} = histogram(n_dat{icond,ilaser}); hold on
        xlim([0 0.5])

        med_f(icond,ilaser) = median(n_dat{icond,ilaser});
        avg_f(icond,ilaser) = mean(n_dat{icond,ilaser});

        h_dat{icond,ilaser}.NumBins = 51;
        h_dat{icond,ilaser}.BinEdges = 0:0.01:0.5;
        h_dat{icond,ilaser}.Normalization = 'probability';
        
        fg = gca;
        ymax = [ymax; fg.YLim(2)];
    end
end
ymax = max(ymax);
cnt = 0;
for icond = 1:8
    for ilaser = 1:2
        cnt = cnt + 1;
        subplot(8,2,cnt);

        ylim([0 ymax])

        hold on
        tmp_v = vline(med_f(icond,ilaser),'r-');
        tmp_v.LineWidth = 1.5;
        hold on
        tmp_v = vline(avg_f(icond,ilaser),'y-');

        if ilaser == 1
            ylabel(cond_list{icond});
        end
        if icond == 1
            if icond == 1
                title('Loff');
            elseif icond == 2
                title('Lon');
            end
        end
        if icond == 8
            xlabel('success/app')
        end
    end
end

if sav_on == 1
    F = getframe(f);
    imwrite(F.cdata, 'ex_success_app_iter1.png', 'png')
    close
end


%%
% for intercepts:
n_dat = [];
for icond = 1:8
    for ilaser = 1:2
        n_dat{icond,ilaser} = [];
    end
end

for icond = 1:8
    for ilaser = 1:2
        n_trials = length(capt_bin{icond,ilaser});
        int_num_tmp = int_num{icond,ilaser};
        capt_bin_tmp = capt_bin{icond,ilaser};

        for iter = 1:n_pts

            try_i = randi(n_trials,1,bstrap_num);

            tot_app_num = sum(int_num_tmp(try_i));
            tot_capt_bin = sum(capt_bin_tmp(try_i));

            % success/num-app
            s_rate = tot_capt_bin/tot_app_num;

            n_dat{icond,ilaser} = [n_dat{icond,ilaser}; s_rate];

            if rem(iter,100) == 0
                if ilaser == 1
                    laser_txt = 'laser off';
                elseif ilaser == 2
                    laser_txt = 'laser on';
                end
                disp(['cond: ' num2str(icond) ', ilaser: ' laser_txt ', num iter: ' num2str(iter) ' of ' num2str(n_pts)])
            end
        end

    end
end
int_dat = n_dat;

g = figure(3202); clf
g.Position = [854 94 984 830]; % 1 monitor setup
cnt = 0;
h_dat = [];
ymax = []; % used to set the ylim afterwards
med_f = [];
avg_f = [];
for icond = 1:8
    for ilaser = 1:2
        cnt = cnt + 1;
        subplot(8,2,cnt); cla
        h_dat{icond,ilaser} = histogram(n_dat{icond,ilaser}); hold on
        xlim([0 1])
        h_dat{icond,ilaser}.NumBins = 41;
        h_dat{icond,ilaser}.BinEdges = 0:.025:1;
        h_dat{icond,ilaser}.Normalization = 'probability';
        
        med_f(icond,ilaser) = median(n_dat{icond,ilaser});
        avg_f(icond,ilaser) = mean(n_dat{icond,ilaser});
        
        fg = gca;
        ymax = [ymax; fg.YLim(2)];

    end
end
ymax = max(ymax);
cnt = 0;
for icond = 1:8
    for ilaser = 1:2
        cnt = cnt + 1;
        subplot(8,2,cnt);

        ylim([0 ymax])

        hold on
        tmp_v = vline(med_f(icond,ilaser),'r-');
        tmp_v.LineWidth = 1.5;
        hold on
        tmp_v = vline(avg_f(icond,ilaser),'y-');

        if ilaser == 1
            ylabel(cond_list{icond});
        end
        if icond == 1
            if icond == 1
                title('Loff');
            elseif icond == 2
                title('Lon');
            end
        end
        if icond == 8
            xlabel('success/int')
        end
    end
end

if sav_on == 1
    F = getframe(g);
    imwrite(F.cdata, 'ex_success_int_iter1.png', 'png')
    close
end

%% do the simulation for 50 experiments each, grouping 5 at a time
% run about 1000 50-exp's? grab the means and medians for each
t_iter = 1000;

bstrap_num = 10; % how many trials
n_pts = 50; % repeat this many

% use 'capt_bin' as it already calculated capture or not (under 1800 fr; laser period)

%%
% for approach:
avg_out = [];
med_out = [];
for icond = 1:8
    for ilaser = 1:2
        avg_out{icond,ilaser} = [];
        med_out{icond,ilaser} = [];
    end
end

for miter = 1:t_iter

    disp([num2str(miter) ' of ' num2str(t_iter)])

    n_dat = [];
    for icond = 1:8
        for ilaser = 1:2
            n_dat{icond,ilaser} = [];
        end
    end
    for icond = 1:8
        for ilaser = 1:2
            n_trials = length(capt_bin{icond,ilaser});
            app_num_tmp = app_num{icond,ilaser};
            capt_bin_tmp = capt_bin{icond,ilaser};
    
            for iter = 1:n_pts
    
                try_i = randi(n_trials,1,bstrap_num);
    
                tot_app_num = sum(app_num_tmp(try_i));
                tot_capt_bin = sum(capt_bin_tmp(try_i));
    
                % success/num-app
                s_rate = tot_capt_bin/tot_app_num;
    
                n_dat{icond,ilaser} = [n_dat{icond,ilaser}; s_rate];
    
            end
    
        end
    end

    for icond = 1:8
        for ilaser = 1:2
            tmp_avg = mean(n_dat{icond,ilaser});
            tmp_med = median(n_dat{icond,ilaser});

            avg_out{icond,ilaser} = [avg_out{icond,ilaser}; tmp_avg];
            med_out{icond,ilaser} = [med_out{icond,ilaser}; tmp_med];
        end
    end

end


% for avg
g = figure(4202); clf
g.Position = [854 94 984 830]; % 1 monitor setup
cnt = 0;
h_dat = [];
ymax = []; % used to set the ylim afterwards
med_f = [];
avg_f = [];
for icond = 1:8
    for ilaser = 1:2
        cnt = cnt + 1;
        subplot(8,2,cnt); cla
        h_dat{icond,ilaser} = histogram(avg_out{icond,ilaser}); hold on
        xlim([0 0.2])
        h_dat{icond,ilaser}.NumBins = 81;
        h_dat{icond,ilaser}.BinEdges = 0:.0025:0.2;
        h_dat{icond,ilaser}.Normalization = 'probability';
        
        fg = gca;
        ymax = [ymax; fg.YLim(2)];

        med_f(icond,ilaser) = median(avg_out{icond,ilaser}); hold on
        avg_f(icond,ilaser) = mean(avg_out{icond,ilaser}); hold on

    end
end
ymax = max(ymax);
cnt = 0;
for icond = 1:8
    for ilaser = 1:2
        cnt = cnt + 1;
        subplot(8,2,cnt);

        ylim([0 ymax])

        hold on
        tmp_v = vline(med_f(icond,ilaser),'r-');
        tmp_v.LineWidth = 1.5;
        hold on
        tmp_v = vline(avg_f(icond,ilaser),'y-');

        if ilaser == 1
            ylabel(cond_list{icond});
        end
        if icond == 1
            if icond == 1
                title('Loff');
            elseif icond == 2
                title('Lon');
            end
        end
        if icond == 8
            xlabel('avg success/app')
        end
    end
end

if sav_on == 1
    F = getframe(g);
    imwrite(F.cdata, '1000iter_success_app_avg.png', 'png')
    close
end


% for median
g = figure(4213); clf
g.Position = [854 94 984 830]; % 1 monitor setup
cnt = 0;
h_dat = [];
ymax = []; % used to set the ylim afterwards
med_f = [];
avg_f = [];
for icond = 1:8
    for ilaser = 1:2
        cnt = cnt + 1;
        subplot(8,2,cnt); cla
        h_dat{icond,ilaser} = histogram(med_out{icond,ilaser}); hold on
        xlim([0 0.2])
        h_dat{icond,ilaser}.NumBins = 81;
        h_dat{icond,ilaser}.BinEdges = 0:.0025:0.2;
        h_dat{icond,ilaser}.Normalization = 'probability';
        
        fg = gca;
        ymax = [ymax; fg.YLim(2)];

        med_f(icond,ilaser) = median(med_out{icond,ilaser}); hold on
        avg_f(icond,ilaser) = mean(med_out{icond,ilaser}); hold on

    end
end
ymax = max(ymax);
cnt = 0;
for icond = 1:8
    for ilaser = 1:2
        cnt = cnt + 1;
        subplot(8,2,cnt);

        ylim([0 ymax])

        hold on
        tmp_v = vline(med_f(icond,ilaser),'r-');
        tmp_v.LineWidth = 1.5;
        hold on
        tmp_v = vline(avg_f(icond,ilaser),'y-');

        if ilaser == 1
            ylabel(cond_list{icond});
        end
        if icond == 1
            if icond == 1
                title('Loff');
            elseif icond == 2
                title('Lon');
            end
        end
        if icond == 8
            xlabel('avg success/app')
        end
    end
end
if sav_on == 1
    F = getframe(g);
    imwrite(F.cdata, '1000iter_success_app_med.png', 'png')
    close
end

%% for intercept:
% do the simulation for 50 experiments each, grouping 5 at a time
% run about 1000 50-exp's? grab the means and medians for each
t_iter = 1000;

bstrap_num = 10; % how many trials
n_pts = 50; % repeat this many

% use 'capt_bin' as it already calculated capture or not (under 1800 fr; laser period)

%%
% for intercept:
avg_out = [];
med_out = [];
for icond = 1:8
    for ilaser = 1:2
        avg_out{icond,ilaser} = [];
        med_out{icond,ilaser} = [];
    end
end

for miter = 1:t_iter

    disp([num2str(miter) ' of ' num2str(t_iter)])

    n_dat = [];
    for icond = 1:8
        for ilaser = 1:2
            n_dat{icond,ilaser} = [];
        end
    end
    for icond = 1:8
        for ilaser = 1:2
            n_trials = length(capt_bin{icond,ilaser});
            int_num_tmp = int_num{icond,ilaser};
            capt_bin_tmp = capt_bin{icond,ilaser};
    
            for iter = 1:n_pts
    
                try_i = randi(n_trials,1,bstrap_num);
    
                tot_app_num = sum(int_num_tmp(try_i));
                tot_capt_bin = sum(capt_bin_tmp(try_i));
    
                % success/num-app
                s_rate = tot_capt_bin/tot_app_num;
    
                n_dat{icond,ilaser} = [n_dat{icond,ilaser}; s_rate];
    
            end
    
        end
    end

    for icond = 1:8
        for ilaser = 1:2
            tmp_avg = mean(n_dat{icond,ilaser},'omitnan');
            tmp_med = median(n_dat{icond,ilaser},'omitnan');

            avg_out{icond,ilaser} = [avg_out{icond,ilaser}; tmp_avg];
            med_out{icond,ilaser} = [med_out{icond,ilaser}; tmp_med];
        end
    end

end


% for avg
g = figure(5212); clf
g.Position = [854 94 984 830]; % 1 monitor setup
cnt = 0;
h_dat = [];
ymax = []; % used to set the ylim afterwards
med_f = [];
avg_f = [];
for icond = 1:8
    for ilaser = 1:2
        cnt = cnt + 1;
        subplot(8,2,cnt); cla
        h_dat{icond,ilaser} = histogram(avg_out{icond,ilaser}); hold on
        xlim([0 0.5])
        h_dat{icond,ilaser}.NumBins = 101;
        h_dat{icond,ilaser}.BinEdges = 0:.005:0.5;
        h_dat{icond,ilaser}.Normalization = 'probability';
        
        fg = gca;
        ymax = [ymax; fg.YLim(2)];

        med_f(icond,ilaser) = median(avg_out{icond,ilaser}); hold on
        avg_f(icond,ilaser) = mean(avg_out{icond,ilaser}); hold on

    end
end
ymax = max(ymax);
cnt = 0;
for icond = 1:8
    for ilaser = 1:2
        cnt = cnt + 1;
        subplot(8,2,cnt);

        ylim([0 ymax])

        hold on
        tmp_v = vline(med_f(icond,ilaser),'r-');
        tmp_v.LineWidth = 1.5;
        hold on
        tmp_v = vline(avg_f(icond,ilaser),'y-');

        if ilaser == 1
            ylabel(cond_list{icond});
        end
        if icond == 1
            if icond == 1
                title('Loff');
            elseif icond == 2
                title('Lon');
            end
        end
        if icond == 8
            xlabel('avg success/int')
        end
    end
end
if sav_on == 1
    F = getframe(g);
    imwrite(F.cdata, '1000iter_success_int_avg.png', 'png')
    close
end

% for median
g = figure(4213); clf
g.Position = [854 94 984 830]; % 1 monitor setup
cnt = 0;
h_dat = [];
ymax = []; % used to set the ylim afterwards
med_f = [];
avg_f = [];
for icond = 1:8
    for ilaser = 1:2
        cnt = cnt + 1;
        subplot(8,2,cnt); cla
        h_dat{icond,ilaser} = histogram(med_out{icond,ilaser}); hold on
        xlim([0 0.5])
        h_dat{icond,ilaser}.NumBins = 201;
        h_dat{icond,ilaser}.BinEdges = 0:.0025:0.5;
        h_dat{icond,ilaser}.Normalization = 'probability';
        
        fg = gca;
        ymax = [ymax; fg.YLim(2)];

        med_f(icond,ilaser) = median(med_out{icond,ilaser}); hold on
        avg_f(icond,ilaser) = mean(med_out{icond,ilaser}); hold on

    end
end
ymax = max(ymax);
cnt = 0;
for icond = 1:8
    for ilaser = 1:2
        cnt = cnt + 1;
        subplot(8,2,cnt);

        ylim([0 ymax])

        hold on
        tmp_v = vline(med_f(icond,ilaser),'r-');
        tmp_v.LineWidth = 1.5;
        hold on
        tmp_v = vline(avg_f(icond,ilaser),'y-');

        if ilaser == 1
            ylabel(cond_list{icond});
        end
        if icond == 1
            if icond == 1
                title('Loff');
            elseif icond == 2
                title('Lon');
            end
        end
        if icond == 8
            xlabel('med success/int')
        end
    end
end

if sav_on == 1
    F = getframe(g);
    imwrite(F.cdata, '1000iter_success_int_med.png', 'png')
    close
end

%% create similar plots for %_of_app_leading_to_int
% use: int_per_app_stat{icond,ilaser} (which has: [num_app num_int bin_capt=(0 or 1)])

t_iter = 1000;
bstrap_num = 10; % how many trials
n_pts = 50; % repeat this many

%
avg_out = [];
med_out = [];
for icond = 1:8
    for ilaser = 1:2
        avg_out{icond,ilaser} = [];
        med_out{icond,ilaser} = [];
    end
end

for miter = 1:t_iter

    disp([num2str(miter) ' of ' num2str(t_iter)])

    n_dat = [];
    for icond = 1:8
        for ilaser = 1:2
            n_dat{icond,ilaser} = [];
        end
    end
    for icond = 1:8
        for ilaser = 1:2
            n_trials = length(capt_bin{icond,ilaser});
            app_num_tmp = int_per_app_stat{icond,ilaser}(:,1);
            int_num_tmp = int_per_app_stat{icond,ilaser}(:,2);
    
            for iter = 1:n_pts
    
                try_i = randi(n_trials,1,bstrap_num);
    
                tot_app_num = sum(app_num_tmp(try_i));
                tot_int_num = sum(int_num_tmp(try_i));
    
                % success/num-app
                s_rate = tot_int_num/tot_app_num;
    
                n_dat{icond,ilaser} = [n_dat{icond,ilaser}; s_rate];
    
            end
    
        end
    end

    for icond = 1:8
        for ilaser = 1:2
            tmp_avg = mean(n_dat{icond,ilaser});
            tmp_med = median(n_dat{icond,ilaser});

            avg_out{icond,ilaser} = [avg_out{icond,ilaser}; tmp_avg];
            med_out{icond,ilaser} = [med_out{icond,ilaser}; tmp_med];
        end
    end

end

% for avg
g = figure(6212); clf
g.Position = [854 94 984 830]; % 1 monitor setup
cnt = 0;
h_dat = [];
ymax = []; % used to set the ylim afterwards
med_f = [];
avg_f = [];
for icond = 1:8
    for ilaser = 1:2
        cnt = cnt + 1;
        subplot(8,2,cnt); cla
        h_dat{icond,ilaser} = histogram(avg_out{icond,ilaser}); hold on
        xlim([0 1])
        h_dat{icond,ilaser}.NumBins = 201;
        h_dat{icond,ilaser}.BinEdges = 0:.005:1;
        h_dat{icond,ilaser}.Normalization = 'probability';
        
        fg = gca;
        ymax = [ymax; fg.YLim(2)];

        med_f(icond,ilaser) = median(avg_out{icond,ilaser}); hold on
        avg_f(icond,ilaser) = mean(avg_out{icond,ilaser}); hold on

    end
end
ymax = max(ymax);
cnt = 0;
for icond = 1:8
    for ilaser = 1:2
        cnt = cnt + 1;
        subplot(8,2,cnt);

        ylim([0 ymax])

        hold on
        tmp_v = vline(med_f(icond,ilaser),'r-');
        tmp_v.LineWidth = 1.5;
        hold on
        tmp_v = vline(avg_f(icond,ilaser),'y-');

        if ilaser == 1
            ylabel(cond_list{icond});
        end
        if icond == 1
            if icond == 1
                title('Loff');
            elseif icond == 2
                title('Lon');
            end
        end
        if icond == 8
            xlabel('avg int/app')
        end
    end
end

if sav_on == 1
    F = getframe(g);
    imwrite(F.cdata, '1000iter_int_per_app_avg.png', 'png')
    close
end

% for median
g = figure(6213); clf
g.Position = [854 94 984 830]; % 1 monitor setup
cnt = 0;
h_dat = [];
ymax = []; % used to set the ylim afterwards
med_f = [];
avg_f = [];
for icond = 1:8
    for ilaser = 1:2
        cnt = cnt + 1;
        subplot(8,2,cnt); cla
        h_dat{icond,ilaser} = histogram(med_out{icond,ilaser}); hold on
        xlim([0 1])
        h_dat{icond,ilaser}.NumBins = 201;
        h_dat{icond,ilaser}.BinEdges = 0:.005:1;
        h_dat{icond,ilaser}.Normalization = 'probability';
        
        fg = gca;
        ymax = [ymax; fg.YLim(2)];

        med_f(icond,ilaser) = median(med_out{icond,ilaser}); hold on
        avg_f(icond,ilaser) = mean(med_out{icond,ilaser}); hold on

    end
end
ymax = max(ymax);
cnt = 0;
for icond = 1:8
    for ilaser = 1:2
        cnt = cnt + 1;
        subplot(8,2,cnt);

        ylim([0 ymax])

        hold on
        tmp_v = vline(med_f(icond,ilaser),'r-');
        tmp_v.LineWidth = 1.5;
        hold on
        tmp_v = vline(avg_f(icond,ilaser),'y-');

        if ilaser == 1
            ylabel(cond_list{icond});
        end
        if icond == 1
            if icond == 1
                title('Loff');
            elseif icond == 2
                title('Lon');
            end
        end
        if icond == 8
            xlabel('med int/app')
        end
    end
end

if sav_on == 1
    F = getframe(g);
    imwrite(F.cdata, '1000iter_int_per_app_med.png', 'png')
    close
end

% example of 1 run
g = figure(6331); clf
g.Position = [854 94 984 830]; % 1 monitor setup
cnt = 0;
h_dat = [];
ymax = []; % used to set the ylim afterwards
med_f = [];
avg_f = [];
for icond = 1:8
    for ilaser = 1:2
        cnt = cnt + 1;
        subplot(8,2,cnt); cla
        h_dat{icond,ilaser} = histogram(n_dat{icond,ilaser}); hold on
        xlim([0 1])
        h_dat{icond,ilaser}.NumBins = 201;
        h_dat{icond,ilaser}.BinEdges = 0:.005:1;
        h_dat{icond,ilaser}.Normalization = 'probability';
        
        med_f(icond,ilaser) = median(n_dat{icond,ilaser});
        avg_f(icond,ilaser) = mean(n_dat{icond,ilaser});
        
        fg = gca;
        ymax = [ymax; fg.YLim(2)];

    end
end
ymax = max(ymax);
cnt = 0;
for icond = 1:8
    for ilaser = 1:2
        cnt = cnt + 1;
        subplot(8,2,cnt);

        ylim([0 ymax])

        hold on
        tmp_v = vline(med_f(icond,ilaser),'r-');
        tmp_v.LineWidth = 1.5;
        hold on
        tmp_v = vline(avg_f(icond,ilaser),'y-');

        if ilaser == 1
            ylabel(cond_list{icond});
        end
        if icond == 1
            if icond == 1
                title('Loff');
            elseif icond == 2
                title('Lon');
            end
        end
        if icond == 8
            xlabel('int/app')
        end
    end
end

if sav_on == 1
    F = getframe(g);
    imwrite(F.cdata, 'ex_int_per_app_iter1.png', 'png')
    close
end