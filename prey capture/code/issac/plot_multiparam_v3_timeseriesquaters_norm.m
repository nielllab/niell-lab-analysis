% 2025_01_05: corrected for app_d_v15 formats

% same as plot_multiparam_v2_timeseriesquaters but this normalizes the
% x-axis into 0 to 1 time-scale

% divide the approaches vs intercepts
% look at from time domain (divide into 4-quaters of min-to-max)
% grab the parameters, put them into the 4-quaters bins

% first, grab the time lengths/distribution

clear all

addpath('f:\prey_capture_lkhd_analysis\code')  
addpath('f:\code')  
addpath('f:\code\hline_vline\')
addpath('f:\code\povilaskarvelis\daboxplot\')
addpath('F:\code\povilaskarvelis\dabarplot\')
addpath('f:\code\densityScatterChart-1.2.0.0\')
addpath('F:\code\shadedErrorBar\')

% load 'app_d'
% load('f:\prey_capture_lkhd_analysis\d_save\app_d_v9.mat')
% load('f:\prey_capture_lkhd_analysis\d_save\app_d_v11.mat')
% load('f:\prey_capture_lkhd_analysis\d_save\app_d_v13.mat')
load('f:\prey_capture_lkhd_analysis\d_save\app_d_v15.mat')
og_app_int = app_int;
app_int = [];
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

% new output:
stack_spd = [];
stack_az = [];
stack_distd = [];
stack_cspd = [];
for icond = 1:8
    for ilaser = 1:2
        stack_spd{icond,ilaser} = [];
        stack_az{icond,ilaser} = [];
        stack_distd{icond,ilaser} = [];
        stack_cspd{icond,ilaser} = [];
        for isuccess = 1:2 % 1=approach(fail), 2=intercept(success)
            stack_spd{icond,ilaser}{isuccess} = [];
            stack_az{icond,ilaser}{isuccess} = [];
            stack_distd{icond,ilaser}{isuccess} = [];
            stack_cspd{icond,ilaser}{isuccess} = [];
        end
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
            s_app = size(tmp_appd,1);
            app_int = [];
            app_int = og_app_int{file_id}(1:s_app);
            % app_int = [];
            % for iapp = 1:size(tmp_appd,1)
            %     st = tmp_appd(iapp,1);
            %     ed = tmp_appd(iapp,2);
            % 
            %     tmp_d2c = p_d2c(st:ed);
            %     f4  = find(tmp_d2c < 4);
            %     if isempty(f4) == 0
            %         app_int = [app_int; 1]; % intercept happened
            %     elseif isempty(f4) == 1
            %         app_int = [app_int; 0]; % intercept did not happen
            %     end
            % end
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

        % create 1800 fr nan array
        tmp_nan_array = nan(1,1800);
        
        if isempty(tmp_appd) == 0 % approaches present
            % use 'app_int' which previously defined (0 = approach, 1 = intercept)
            t_app = size(tmp_appd,1);
            t_app_int = length(app_int);
            if t_app ~= t_app_int
                disp(['Something wrong with: #' num2str(ifile) ' [does not match: the appd and app_int numbers]'])
            end
            for iapp = 1:t_app
                
                st = tmp_appd(iapp,1);
                ed = tmp_appd(iapp,2);

                n_nan_spd = tmp_nan_array;
                n_nan_az = tmp_nan_array;
                n_nan_distd = tmp_nan_array;
                n_nan_cspd = tmp_nan_array;

                tmp_spd = app_p_nan{ifile}(1,st:ed);
                % tmp_spd = app_p{ifile}(1,st:ed);
                sav_len = length(tmp_spd);
                n_nan_spd(1:sav_len) = tmp_spd;

                tmp_az = app_p_nan{ifile}(2,st:ed);
                % tmp_az = app_p{ifile}(2,st:ed);
                n_nan_az(1:sav_len) = tmp_az;

                tmp_distd = app_p_nan{ifile}(3,st:ed);
                % tmp_distd = app_p{ifile}(3,st:ed);
                n_nan_distd(1:sav_len) = tmp_distd;

                % tmp_cspd = app_p_nan{ifile}(4,st:ed);
                tmp_cspd = app_p{ifile}(4,st:ed);
                n_nan_cspd(1:sav_len) = tmp_cspd;

                % save into the 'stack_spd', 'stack_az', 'stack_distd',
                % 'stack_cspd'
                int_loc = app_int(iapp);
                if int_loc == 0
                    stack_spd{cond_loc,laser_loc}{1} = [stack_spd{cond_loc,laser_loc}{1}; n_nan_spd];
                    stack_az{cond_loc,laser_loc}{1} = [stack_az{cond_loc,laser_loc}{1}; n_nan_az];
                    stack_distd{cond_loc,laser_loc}{1} = [stack_distd{cond_loc,laser_loc}{1}; n_nan_distd];
                    stack_cspd{cond_loc,laser_loc}{1} = [stack_cspd{cond_loc,laser_loc}{1}; n_nan_cspd];
                elseif int_loc == 1
                    stack_spd{cond_loc,laser_loc}{2} = [stack_spd{cond_loc,laser_loc}{2}; n_nan_spd];
                    stack_az{cond_loc,laser_loc}{2} = [stack_az{cond_loc,laser_loc}{2}; n_nan_az];
                    stack_distd{cond_loc,laser_loc}{2} = [stack_distd{cond_loc,laser_loc}{2}; n_nan_distd];
                    stack_cspd{cond_loc,laser_loc}{2} = [stack_cspd{cond_loc,laser_loc}{2}; n_nan_cspd];
                end


            end

        % elseif isempty(tmp_appd) == 1 % approaches absent
            % do nothing
        end

    end
    % end
end


% - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - %
% - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - %
% - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - %
% - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - %
% - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - %
% - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - %
% - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - %
% - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - %
% - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - %
% - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - %
% - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - % - %

% Organize and plot
% p1_spd = [];
% p1_az = [];
% p1_dist = [];
% p1_cspd = [];

% cd('F:\prey_capture_lkhd_analysis\figure\2024_12_18\')
cd('F:\prey_capture_lkhd_analysis\figure\2025_01_06b\')

% what is the max y-axis so you can normalize across all plots
max_spd = [];
max_az = [];
max_dist = [];
max_cspd = [];

for icond = 1:8
    for ilaser = 1:2
        % p1_spd{icond,ilaser} = [];
        % p1_az{icond,ilaser} = [];
        % p1_dist{icond,ilaser} = [];
        % p1_cspd{icond,ilaser} = [];
        for isuccess = 1:2
            % p1_spd{icond,ilaser}{isuccess} = [];
            % p1_az{icond,ilaser}{isuccess} = [];
            % p1_dist{icond,ilaser}{isuccess} = [];
            % p1_cspd{icond,ilaser}{isuccess} = [];
            
            max_spd = [max_spd; (max(stack_spd{icond,ilaser}{isuccess}))];
            max_az = [max_az; (max(stack_az{icond,ilaser}{isuccess}))];
            max_dist = [max_dist; (max(stack_distd{icond,ilaser}{isuccess}))];
            max_cspd = [max_cspd; (max(stack_cspd{icond,ilaser}{isuccess}))];
        end

        
    end
end

pmax_spd = [];
pmax_az = [];
pmax_dist = [];
pmax_cspd = [];

for icond = 1:8
    for ilaser = 1:2
        for isuccess = 1:2

            fisnan = isnan(stack_spd{icond,ilaser}{isuccess}(:));
            tmp = stack_spd{icond,ilaser}{isuccess}(:);
            tmp(fisnan) = [];
            pmax_spd = [pmax_spd; tmp];

            fisnan = isnan(stack_az{icond,ilaser}{isuccess}(:));
            tmp = stack_az{icond,ilaser}{isuccess}(:);
            tmp(fisnan) = [];
            pmax_az = [pmax_az; tmp];

            fisnan = isnan(stack_distd{icond,ilaser}{isuccess}(:));
            tmp = stack_distd{icond,ilaser}{isuccess}(:);
            tmp(fisnan) = [];
            pmax_dist = [pmax_dist; tmp];

            fisnan = isnan(stack_cspd{icond,ilaser}{isuccess}(:));
            tmp = stack_cspd{icond,ilaser}{isuccess}(:);
            tmp(fisnan) = [];
            pmax_cspd = [pmax_cspd; tmp];
        end
    end
end

% calculate the mean and 2 std's so you can correctly plot the ranges
pstat_spd = [mean(pmax_spd) std(pmax_spd) mean(pmax_spd)+2*std(pmax_spd)];
pstat_az = [mean(pmax_az) std(pmax_az) mean(pmax_az)+2*std(pmax_az)];
pstat_dist = [mean(pmax_dist) std(pmax_dist) mean(pmax_dist)+2*std(pmax_dist)];
pstat_cspd = [mean(pmax_cspd) std(pmax_cspd) mean(pmax_cspd)+2*std(pmax_cspd)];

% spd
imp = [];
for icond = 1:8
    f = figure(2000+icond); clf
    f.Position = [658 148 1205 745];
    ylimmax = [];
    xlimmax = [];
    climmax = [];

    stack_p = [];
    for ilaser = 1:2
        imp{icond,ilaser} = [];
        for isuccess = 1:2
            subplot(2,3,(ilaser-1)*3+isuccess), cla

            stack_p = [stack_p; stack_spd{icond,ilaser}{isuccess}];

            tmp = mean(stack_spd{icond,ilaser}{isuccess},'omitnan');
            fnan = isnan(tmp);
            f1 = find(fnan == 1);
            f1 = f1(1);

            tmp = stack_spd{icond,ilaser}{isuccess}(:,1:f1);

            imp{icond,ilaser}{isuccess} = imagesc(tmp);

            ax1 = gca;
            ylimmax = [ylimmax; ax1.YLim(2)];
            xlimmax = [xlimmax; ax1.XLim(2)];
            climmax = [climmax; ax1.CLim(2)];

            if ilaser == 1
                lasertext = 'Loff';
            elseif ilaser == 2
                lasertext = 'Lon';
            end
            if isuccess == 1
                successtext = 'approach';
            elseif isuccess == 2
                successtext = 'intercept';
            end

            t_text = [cond_list{icond} ', ' successtext ', ' lasertext];
            title(t_text)

            colorbar
        end
    end
    for ilaser = 1:2
        for isuccess = 1:2
            subplot(2,3,(ilaser-1)*3+isuccess)

            ylim([0 max(ylimmax)])
            clim([0 max(climmax)])
        end
    end
    % stacked
    subplot(2,3,[3 6])
    imagesc(stack_p)
    xlim([0 max(xlimmax)])
    colorbar


    F = getframe(f);
    imwrite(F.cdata, [cond_list{icond} '_spd_stacked_cmap.png'], 'png')
    close
end

close all


% line chart:

imp = [];
% cmap = [0.75 0 0; 0 0 0];
cmap = [0 0.75 0.75; 0 0 0];
xx = [];
yy = [];
for icond = 1:8
    f = figure(2100+icond); clf
    f.Position = [1176 85 602 846];
    ylimmax = [];
    xlimmax = [];

    for ilaser = 1:2

        subplot(2,1,ilaser)
        imp{icond,ilaser} = [];
        for isuccess = 1:2

            tmp_p = stack_spd{icond,ilaser}{isuccess};
            tmp = [];
            for itmp = 1:size(tmp_p,1)
                smp = tmp_p(itmp,:);
                
                if sum(isnan(smp)) == 1800
                    % do nothing
                else
                    fnan = isnan(smp);
                    if sum(fnan) < 1800-5
                        df_fnan = diff(fnan);
                        fdf_fnan = find(df_fnan == 1);
                        if length(fdf_fnan) == 1
                            final_pt = fdf_fnan;
                        else
                            final_pt = fdf_fnan(end);
                        end
                        smp = tmp_p(itmp,1:final_pt);
                        len_smp = length(smp);
                        x_smp = 1:((len_smp-1)/(1800-1)):(len_smp);
                        % x_smp = 1:(len_smp/(1800-1)):(len_smp+1);
                        y_interp = interp1(1:len_smp,smp,x_smp);
        
                        tmp = [tmp; y_interp];
                    end
                end
                
            end

            % tmp = mean(stack_spd{icond,ilaser}{isuccess},'omitnan');
            tmp2 = plot(tmp'); hold on
            for itmp = 1:length(tmp2)
                % tmp(itmp).Color = [cmap(isuccess,:) 0.15];
                tmp2(itmp).Color = [cmap(isuccess,:) 0.2];
            end

            ax1 = gca;
            ylimmax = [ylimmax; ax1.YLim(2)];
            xlimmax = [xlimmax; ax1.XLim(2)];

            if ilaser == 1
                lasertext = 'Loff';
            elseif ilaser == 2
                lasertext = 'Lon';
            end

            % if isuccess == 2
            %     legend('fail','success')
            % end

            t_text = [cond_list{icond} ', ' lasertext];
            title(t_text)
        end
    end
    for ilaser = 1:2

        subplot(2,1,ilaser)
        for isuccess = 1:2

            % ylim([0 max(ylimmax)])
            ylim([0 50])
            xlim([0 max(xlimmax)])
        end
    end
    xx = [xx; xlimmax];
    yy = [yy; ylimmax];
    % 
    F = getframe(f);
    imwrite(F.cdata, [cond_list{icond} '_spd_indiv_line_norm.png'], 'png')
    close
end


% % just plot 1 by 1, and don't interpolate
% % 
% imp = [];
% % cmap = [0.75 0 0; 0 0 0];
% cmap = [0 0.75 0.75; 0 0 0];
% xx = [];
% yy = [];
% for icond = 1:8
%     f = figure(2100+icond); clf
%     f.Position = [1176 85 602 846];
%     ylimmax = [];
%     xlimmax = [];
% 
%     for ilaser = 1:2
% 
%         subplot(2,1,ilaser)
%         imp{icond,ilaser} = [];
%         for isuccess = 1:2
% 
%             tmp_p = stack_spd{icond,ilaser}{isuccess};
%             tmp = [];
%             cnt_tmp = 0;
%             for itmp = 1:size(tmp_p,1)
%                 smp = tmp_p(itmp,:);
% 
%                 if sum(isnan(smp)) == 1800
%                     % do nothing
%                 else
%                     fnan = isnan(smp);
%                     if sum(fnan) < 1800-5
%                         df_fnan = diff(fnan);
%                         fdf_fnan = find(df_fnan == 1);
%                         if length(fdf_fnan) == 1
%                             final_pt = fdf_fnan;
%                         else
%                             final_pt = fdf_fnan(end);
%                         end
%                         smp = tmp_p(itmp,1:final_pt);
%                         len_smp = length(smp);
%                         % x_smp = 0:(len_smp/(1800-1)):len_smp;
%                         % y_interp = interp1(1:len_smp,smp,x_smp);
%                         x_smp = 0:(1/(len_smp-1)):1;
% 
%                         cnt_tmp = cnt_tmp + 1;
%                         tmp{cnt_tmp} = plot(x_smp,smp); hold on
% 
%                         % tmp = [tmp; y_interp];
%                     end
%                 end
% 
%             end
% 
%             % tmp = mean(stack_spd{icond,ilaser}{isuccess},'omitnan');
%             % tmp2 = plot(tmp'); hold on
%             % for itmp = 1:length(tmp2)
%             for itmp = 1:length(tmp)
%                 tmp{itmp}.Color = [cmap(isuccess,:) 0.2];
%                 % tmp2(itmp).Color = [cmap(isuccess,:) 0.2];
%             end
% 
%             ax1 = gca;
%             ylimmax = [ylimmax; ax1.YLim(2)];
%             % xlimmax = [xlimmax; ax1.XLim(2)];
% 
%             if ilaser == 1
%                 lasertext = 'Loff';
%             elseif ilaser == 2
%                 lasertext = 'Lon';
%             end
% 
%             % if isuccess == 2
%             %     legend('fail','success')
%             % end
% 
%             t_text = [cond_list{icond} ', ' lasertext];
%             title(t_text)
%         end
%     end
%     for ilaser = 1:2
% 
%         subplot(2,1,ilaser)
%         for isuccess = 1:2
% 
%             % ylim([0 max(ylimmax)])
%             ylim([0 50])
%             % xlim([0 max(xlimmax)])
%         end
%     end
%     % xx = [xx; xlimmax];
%     yy = [yy; ylimmax];
% 
%     % F = getframe(f);
%     % imwrite(F.cdata, [cond_list{icond} '_spd_indiv_line_norm_1by1.png'], 'png')
%     % close
% end
% 
% % average by binning
% % 
% imp = [];
% x_bin = -0.025:0.05:1.025;
% % cmap = [0.75 0 0; 0 0 0];
% cmap = [0 0.75 0.75; 0 0 0];
% xx = [];
% yy = [];
% for icond = 1:8
%     f = figure(2100+icond); clf
%     f.Position = [1176 85 602 846];
%     ylimmax = [];
%     xlimmax = [];
% 
%     for ilaser = 1:2
% 
%         subplot(2,1,ilaser)
%         imp{icond,ilaser} = [];
% 
%         for isuccess = 1:2
% 
%             tmp_p = stack_spd{icond,ilaser}{isuccess};
%             tmp = [];
%             cnt_tmp = 0;
%             for itmp = 1:size(tmp_p,1)
%                 smp = tmp_p(itmp,:);
% 
%                 if sum(isnan(smp)) == 1800
%                     % do nothing
%                 else
%                     fnan = isnan(smp);
%                     if sum(fnan) < 1800-5
%                         df_fnan = diff(fnan);
%                         fdf_fnan = find(df_fnan == 1);
%                         if length(fdf_fnan) == 1
%                             final_pt = fdf_fnan;
%                         else
%                             final_pt = fdf_fnan(end);
%                         end
%                         smp = tmp_p(itmp,1:final_pt);
%                         len_smp = length(smp);
%                         % x_smp = 0:(len_smp/(1800-1)):len_smp;
%                         % y_interp = interp1(1:len_smp,smp,x_smp);
%                         x_smp = 0:(1/(len_smp-1)):1;
% 
%                         % cnt_tmp = cnt_tmp + 1;
%                         % tmp{cnt_tmp} = plot(x_smp,smp); hold on
% 
%                         % tmp = [tmp; y_interp];
% 
%                         for ismp = 1:length(x_smp)
%                             t_smp = x_smp(ismp);
%                             y_smp = smp(ismp);
%                             for ixbin = 1:(length(x_bin)-1)
%                                 xst = x_bin(ixbin);
%                                 xed = x_bin(ixbin+1);
% 
%                                 if xst < t_smp && xed > t_smp
%                                     xloc = find(x_bin == xst);
%                                     imp{icond,ilaser}(itmp,xloc) = y_smp;
%                                 end
%                             end
%                         end
%                     end
%                 end
% 
%             end
% 
%             % tmp = mean(stack_spd{icond,ilaser}{isuccess},'omitnan');
%             % tmp2 = plot(tmp'); hold on
%             tmp2 = plot(imp{icond,ilaser}'); hold on
%             % for itmp = 1:length(tmp2)
%             for itmp = 1:length(tmp)
%                 % tmp{itmp}.Color = [cmap(isuccess,:) 0.2];
%                 tmp2(itmp).Color = [cmap(isuccess,:) 0.2];
%             end
% 
%             ax1 = gca;
%             ylimmax = [ylimmax; ax1.YLim(2)];
%             % xlimmax = [xlimmax; ax1.XLim(2)];
% 
%             if ilaser == 1
%                 lasertext = 'Loff';
%             elseif ilaser == 2
%                 lasertext = 'Lon';
%             end
% 
%             % if isuccess == 2
%             %     legend('fail','success')
%             % end
% 
%             t_text = [cond_list{icond} ', ' lasertext];
%             title(t_text)
%         end
%     end
%     for ilaser = 1:2
% 
%         subplot(2,1,ilaser)
%         for isuccess = 1:2
% 
%             % ylim([0 max(ylimmax)])
%             ylim([0 50])
%             % xlim([0 max(xlimmax)])
%         end
%     end
%     % xx = [xx; xlimmax];
%     yy = [yy; ylimmax];
% 
%     F = getframe(f);
%     imwrite(F.cdata, [cond_list{icond} '_spd_indiv_line_norm_1by1.png'], 'png')
%     close
% end


% avg
imp = [];
% cmap = [0.75 0 0; 0 0 0];
cmap = [0 0.75 0.75; 0 0 0];

xx = [];
yy = [];
for icond = 1:8
    f = figure(2100+icond); clf
    f.Position = [1176 85 602 846];
    ylimmax = [];
    xlimmax = [];

    for ilaser = 1:2

        subplot(2,1,ilaser)
        imp{icond,ilaser} = [];

        tmp = [];
        for isuccess = 1:2

            tmp_p = stack_spd{icond,ilaser}{isuccess};
            tmp = [];
            for itmp = 1:size(tmp_p,1)
                smp = tmp_p(itmp,:);

                if sum(isnan(smp)) == 1800
                    % do nothing
                else
                    fnan = isnan(smp);
                    if sum(fnan) < 1800-5
                        df_fnan = diff(fnan);
                        fdf_fnan = find(df_fnan == 1);
                        if length(fdf_fnan) == 1
                            final_pt = fdf_fnan;
                        else
                            final_pt = fdf_fnan(end);
                        end
                        smp = tmp_p(itmp,1:final_pt);
                        len_smp = length(smp);
                        x_smp = 1:((len_smp-1)/(1800-1)):(len_smp);
                        % x_smp = 1:(len_smp/(1800-1)):(len_smp+1);
                        y_interp = interp1(1:len_smp,smp,x_smp);

                        tmp = [tmp; y_interp];
                    end
                end

            end

            % tmp_avg = mean(stack_spd{icond,ilaser}{isuccess},'omitnan');
            tmp_avg = mean(tmp,'omitnan');
            % tmp_std = [];
            % for icol = 1:size(stack_spd{icond,ilaser}{isuccess},2)
            %     tmp_col = stack_spd{icond,ilaser}{isuccess}(:,icol);
            %     n_nan = isnan(tmp_col);
            %     n_nan = sum(n_nan);
            %     n_num = size(stack_spd{icond,ilaser}{isuccess},1)-n_nan;
            %     tmp_std(icol) = std(tmp_col,'omitnan')./sqrt(n_num);
            % end
            tmp_std = std(tmp,'omitnan')./sqrt(size(tmp,1));

            x_plot = 1:1800;
            y_plot = tmp_avg;
            e_plot = tmp_std;

            s_tmp = shadedErrorBar(x_plot,y_plot,e_plot); hold on
            s_tmp.mainLine.Color = cmap(isuccess,:);
            s_tmp.edge(1).LineStyle = 'none';
            s_tmp.edge(2).LineStyle = 'none';
            s_tmp.patch.FaceColor = cmap(isuccess,:);

            ax1 = gca;
            ylimmax = [ylimmax; ax1.YLim(2)];
            xlimmax = [xlimmax; ax1.XLim(2)];

            if ilaser == 1
                lasertext = 'Loff';
            elseif ilaser == 2
                lasertext = 'Lon';
            end

            % if isuccess == 2
            %     legend('fail','success')
            % end

            t_text = [cond_list{icond} ', ' lasertext];
            title(t_text)
        end
    end
    for ilaser = 1:2

        subplot(2,1,ilaser)
        for isuccess = 1:2

            % ylim([0 max(ylimmax)])
            ylim([0 30])
            xlim([0 max(xlimmax)])
        end
    end

    xx = [xx; xlimmax];
    yy = [yy; ylimmax];

    F = getframe(f);
    imwrite(F.cdata, [cond_list{icond} '_spd_avg_line_norm.png'], 'png')
    close
end














% az
imp = [];
for icond = 1:8
    f = figure(3000+icond); clf
    f.Position = [658 148 1205 745];
    ylimmax = [];
    xlimmax = [];
    climmax = [];

    stack_p = [];
    for ilaser = 1:2
        imp{icond,ilaser} = [];
        for isuccess = 1:2
            subplot(2,3,(ilaser-1)*3+isuccess), cla

            stack_p = [stack_p; stack_az{icond,ilaser}{isuccess}];

            tmp = mean(stack_az{icond,ilaser}{isuccess},'omitnan');
            fnan = isnan(tmp);
            f1 = find(fnan == 1);
            f1 = f1(1);

            tmp = stack_az{icond,ilaser}{isuccess}(:,1:f1);

            imp{icond,ilaser}{isuccess} = imagesc(tmp);

            ax1 = gca;
            ylimmax = [ylimmax; ax1.YLim(2)];
            xlimmax = [xlimmax; ax1.XLim(2)];
            climmax = [climmax; ax1.CLim(2)];

            if ilaser == 1
                lasertext = 'Loff';
            elseif ilaser == 2
                lasertext = 'Lon';
            end
            if isuccess == 1
                successtext = 'approach';
            elseif isuccess == 2
                successtext = 'intercept';
            end

            t_text = [cond_list{icond} ', ' successtext ', ' lasertext];
            title(t_text)

            colorbar
        end
    end
    for ilaser = 1:2
        for isuccess = 1:2
            subplot(2,3,(ilaser-1)*3+isuccess)

            ylim([0 max(ylimmax)])
            clim([0 max(climmax)])
        end
    end
    % stacked
    subplot(2,3,[3 6])
    imagesc(stack_p)
    xlim([0 max(xlimmax)])
    colorbar


    F = getframe(f);
    imwrite(F.cdata, [cond_list{icond} '_az_stacked_cmap.png'], 'png')
    close
end

% line chart:

imp = [];
xx = [];
yy = [];
% cmap = [0.75 0 0; 0 0 0];
cmap = [0 0.75 0.75; 0 0 0];
for icond = 1:8
    f = figure(3100+icond); clf
    f.Position = [1176 85 602 846];
    ylimmax = [];
    xlimmax = [];

    for ilaser = 1:2

        subplot(2,1,ilaser)
        imp{icond,ilaser} = [];
        for isuccess = 1:2

            tmp_p = stack_az{icond,ilaser}{isuccess};
            tmp = [];
            for itmp = 1:size(tmp_p,1)
                smp = tmp_p(itmp,:);
                
                if sum(isnan(smp)) == 1800
                    % do nothing
                else
                    fnan = isnan(smp);
                    if sum(fnan) < 1800-5
                        df_fnan = diff(fnan);
                        fdf_fnan = find(df_fnan == 1);
                        if length(fdf_fnan) == 1
                            final_pt = fdf_fnan;
                        else
                            final_pt = fdf_fnan(end);
                        end
                        smp = tmp_p(itmp,1:final_pt);
                        len_smp = length(smp);
                        x_smp = 1:((len_smp-1)/(1800-1)):(len_smp);
                        % x_smp = 1:(len_smp/(1800-1)):(len_smp+1);
                        y_interp = interp1(1:len_smp,smp,x_smp);
        
                        tmp = [tmp; y_interp];
                    end
                end
                
            end

            % tmp = mean(stack_az{icond,ilaser}{isuccess},'omitnan');
            tmp2 = plot(tmp'); hold on
            for itmp = 1:length(tmp2)
                % tmp(itmp).Color = [cmap(isuccess,:) 0.15];
                tmp2(itmp).Color = [cmap(isuccess,:) 0.2];
            end

            ax1 = gca;
            ylimmax = [ylimmax; ax1.YLim(2)];
            xlimmax = [xlimmax; ax1.XLim(2)];

            if ilaser == 1
                lasertext = 'Loff';
            elseif ilaser == 2
                lasertext = 'Lon';
            end

            % if isuccess == 2
            %     legend('fail','success')
            % end

            t_text = [cond_list{icond} ', ' lasertext];
            title(t_text)
        end
    end
    for ilaser = 1:2

        subplot(2,1,ilaser)
        for isuccess = 1:2

            % ylim([0 max(ylimmax)])
            ylim([0 60])
            xlim([0 max(xlimmax)])
        end
    end
    xx = [xx; xlimmax];
    yy = [yy; ylimmax];

    F = getframe(f);
    imwrite(F.cdata, [cond_list{icond} '_az_indiv_line_norm.png'], 'png')
    close
end


% % just plot 1 by 1, and don't interpolate
% % 
% imp = [];
% % cmap = [0.75 0 0; 0 0 0];
% cmap = [0 0.75 0.75; 0 0 0];
% xx = [];
% yy = [];
% for icond = 1:8
%     f = figure(2100+icond); clf
%     f.Position = [1176 85 602 846];
%     ylimmax = [];
%     xlimmax = [];
% 
%     for ilaser = 1:2
% 
%         subplot(2,1,ilaser)
%         imp{icond,ilaser} = [];
%         for isuccess = 1:2
% 
%             tmp_p = stack_az{icond,ilaser}{isuccess};
%             tmp = [];
%             cnt_tmp = 0;
%             for itmp = 1:size(tmp_p,1)
%                 smp = tmp_p(itmp,:);
% 
%                 if sum(isnan(smp)) == 1800
%                     % do nothing
%                 else
%                     fnan = isnan(smp);
%                     if sum(fnan) < 1800-5
%                         df_fnan = diff(fnan);
%                         fdf_fnan = find(df_fnan == 1);
%                         if length(fdf_fnan) == 1
%                             final_pt = fdf_fnan;
%                         else
%                             final_pt = fdf_fnan(end);
%                         end
%                         smp = tmp_p(itmp,1:final_pt);
%                         len_smp = length(smp);
%                         % x_smp = 0:(len_smp/(1800-1)):len_smp;
%                         % y_interp = interp1(1:len_smp,smp,x_smp);
%                         x_smp = 0:(1/(len_smp-1)):1;
% 
%                         cnt_tmp = cnt_tmp + 1;
%                         tmp{cnt_tmp} = plot(x_smp,smp); hold on
% 
%                         % tmp = [tmp; y_interp];
%                     end
%                 end
% 
%             end
% 
%             % tmp = mean(stack_spd{icond,ilaser}{isuccess},'omitnan');
%             % tmp2 = plot(tmp'); hold on
%             % for itmp = 1:length(tmp2)
%             for itmp = 1:length(tmp)
%                 tmp{itmp}.Color = [cmap(isuccess,:) 0.2];
%                 % tmp2(itmp).Color = [cmap(isuccess,:) 0.2];
%             end
% 
%             ax1 = gca;
%             ylimmax = [ylimmax; ax1.YLim(2)];
%             % xlimmax = [xlimmax; ax1.XLim(2)];
% 
%             if ilaser == 1
%                 lasertext = 'Loff';
%             elseif ilaser == 2
%                 lasertext = 'Lon';
%             end
% 
%             % if isuccess == 2
%             %     legend('fail','success')
%             % end
% 
%             t_text = [cond_list{icond} ', ' lasertext];
%             title(t_text)
%         end
%     end
%     for ilaser = 1:2
% 
%         subplot(2,1,ilaser)
%         for isuccess = 1:2
% 
%             % ylim([0 max(ylimmax)])
%             ylim([0 50])
%             % xlim([0 max(xlimmax)])
%         end
%     end
%     % xx = [xx; xlimmax];
%     yy = [yy; ylimmax];
% 
%     F = getframe(f);
%     imwrite(F.cdata, [cond_list{icond} '_az_indiv_line_norm_1by1.png'], 'png')
%     close
% end

% avg
imp = [];
% cmap = [0.75 0 0; 0 0 0];
cmap = [0 0.75 0.75; 0 0 0];
for icond = 1:8
    f = figure(2100+icond); clf
    f.Position = [1176 85 602 846];
    ylimmax = [];
    xlimmax = [];

    for ilaser = 1:2

        subplot(2,1,ilaser)
        imp{icond,ilaser} = [];

        tmp = [];
        for isuccess = 1:2

            tmp_p = stack_az{icond,ilaser}{isuccess};
            tmp = [];
            for itmp = 1:size(tmp_p,1)
                smp = tmp_p(itmp,:);

                if sum(isnan(smp)) == 1800
                    % do nothing
                else
                    fnan = isnan(smp);
                    if sum(fnan) < 1800-5
                        df_fnan = diff(fnan);
                        fdf_fnan = find(df_fnan == 1);
                        if length(fdf_fnan) == 1
                            final_pt = fdf_fnan;
                        else
                            final_pt = fdf_fnan(end);
                        end
                        smp = tmp_p(itmp,1:final_pt);
                        len_smp = length(smp);
                        x_smp = 1:((len_smp-1)/(1800-1)):(len_smp);
                        % x_smp = 1:(len_smp/(1800-1)):(len_smp+1);
                        y_interp = interp1(1:len_smp,smp,x_smp);

                        tmp = [tmp; y_interp];
                    end
                end

            end

            % tmp_avg = mean(stack_az{icond,ilaser}{isuccess},'omitnan');
            tmp_avg = mean(tmp,'omitnan');
            % tmp_std = [];
            % for icol = 1:size(stack_az{icond,ilaser}{isuccess},2)
            %     tmp_col = stack_az{icond,ilaser}{isuccess}(:,icol);
            %     n_nan = isnan(tmp_col);
            %     n_nan = sum(n_nan);
            %     n_num = size(stack_az{icond,ilaser}{isuccess},1)-n_nan;
            %     tmp_std(icol) = std(tmp_col,'omitnan')./sqrt(n_num);
            % end
            tmp_std = std(tmp,'omitnan')./sqrt(size(tmp,1));

            x_plot = 1:1800;
            y_plot = tmp_avg;
            e_plot = tmp_std;

            s_tmp = shadedErrorBar(x_plot,y_plot,e_plot); hold on
            s_tmp.mainLine.Color = cmap(isuccess,:);
            s_tmp.edge(1).LineStyle = 'none';
            s_tmp.edge(2).LineStyle = 'none';
            s_tmp.patch.FaceColor = cmap(isuccess,:);

            ax1 = gca;
            ylimmax = [ylimmax; ax1.YLim(2)];
            xlimmax = [xlimmax; ax1.XLim(2)];

            if ilaser == 1
                lasertext = 'Loff';
            elseif ilaser == 2
                lasertext = 'Lon';
            end

            % if isuccess == 2
            %     legend('fail','success')
            % end

            t_text = [cond_list{icond} ', ' lasertext];
            title(t_text)
        end
    end
    for ilaser = 1:2

        subplot(2,1,ilaser)
        for isuccess = 1:2

            % ylim([0 max(ylimmax)])
            ylim([0 60])
            xlim([0 max(xlimmax)])
        end
    end

    F = getframe(f);
    imwrite(F.cdata, [cond_list{icond} '_az_avg_line_norm.png'], 'png')
    close
end











% fix the 1by1 part from hereon:
% (done, 2025_01_06)


% distd
imp = [];
for icond = 1:8
    f = figure(4000+icond); clf
    f.Position = [658 148 1205 745];
    ylimmax = [];
    xlimmax = [];
    climmax = [];

    stack_p = [];
    for ilaser = 1:2
        imp{icond,ilaser} = [];
        for isuccess = 1:2
            subplot(2,3,(ilaser-1)*3+isuccess), cla

            stack_p = [stack_p; stack_distd{icond,ilaser}{isuccess}];

            tmp = mean(stack_distd{icond,ilaser}{isuccess},'omitnan');
            fnan = isnan(tmp);
            f1 = find(fnan == 1);
            f1 = f1(1);

            tmp = stack_distd{icond,ilaser}{isuccess}(:,1:f1);

            imp{icond,ilaser}{isuccess} = imagesc(tmp);

            ax1 = gca;
            ylimmax = [ylimmax; ax1.YLim(2)];
            xlimmax = [xlimmax; ax1.XLim(2)];
            climmax = [climmax; ax1.CLim(2)];

            if ilaser == 1
                lasertext = 'Loff';
            elseif ilaser == 2
                lasertext = 'Lon';
            end
            if isuccess == 1
                successtext = 'approach';
            elseif isuccess == 2
                successtext = 'intercept';
            end

            t_text = [cond_list{icond} ', ' successtext ', ' lasertext];
            title(t_text)

            colorbar
        end
    end
    for ilaser = 1:2
        for isuccess = 1:2
            subplot(2,3,(ilaser-1)*3+isuccess)

            ylim([0 max(ylimmax)])
            clim([0 max(climmax)])
        end
    end
    % stacked
    subplot(2,3,[3 6])
    imagesc(stack_p)
    xlim([0 max(xlimmax)])
    colorbar

    F = getframe(f);
    imwrite(F.cdata, [cond_list{icond} '_distd_stacked_cmap.png'], 'png')
    close
end

% line chart:

imp = [];
xx = [];
yy = [];
% cmap = [0.75 0 0; 0 0 0];
cmap = [0 0.75 0.75; 0 0 0];
for icond = 1:8
    f = figure(3100+icond); clf
    f.Position = [1176 85 602 846];
    ylimmax = [];
    xlimmax = [];

    for ilaser = 1:2

        subplot(2,1,ilaser)
        imp{icond,ilaser} = [];
        for isuccess = 1:2

            tmp_p = stack_distd{icond,ilaser}{isuccess};
            tmp = [];
            for itmp = 1:size(tmp_p,1)
                smp = tmp_p(itmp,:);
                
                if sum(isnan(smp)) == 1800
                    % do nothing
                else
                    fnan = isnan(smp);
                    if sum(fnan) < 1800-5
                        df_fnan = diff(fnan);
                        fdf_fnan = find(df_fnan == 1);
                        if length(fdf_fnan) == 1
                            final_pt = fdf_fnan;
                        else
                            final_pt = fdf_fnan(end);
                        end
                        smp = tmp_p(itmp,1:final_pt);
                        len_smp = length(smp);
                        x_smp = 1:((len_smp-1)/(1800-1)):(len_smp);
                        % x_smp = 1:(len_smp/(1800-1)):(len_smp+1);
                        y_interp = interp1(1:len_smp,smp,x_smp);
        
                        tmp = [tmp; y_interp];
                    end
                end
                
            end

            % tmp = mean(stack_distd{icond,ilaser}{isuccess},'omitnan');
            tmp2 = plot(tmp'); hold on
            for itmp = 1:length(tmp2)
                % tmp(itmp).Color = [cmap(isuccess,:) 0.15];
                tmp2(itmp).Color = [cmap(isuccess,:) 0.2];
            end

            ax1 = gca;
            ylimmax = [ylimmax; ax1.YLim(2)];
            xlimmax = [xlimmax; ax1.XLim(2)];

            if ilaser == 1
                lasertext = 'Loff';
            elseif ilaser == 2
                lasertext = 'Lon';
            end

            % if isuccess == 2
            %     legend('fail','success')
            % end

            t_text = [cond_list{icond} ', ' lasertext];
            title(t_text)
        end
    end
    for ilaser = 1:2

        subplot(2,1,ilaser)
        for isuccess = 1:2

            % ylim([0 max(ylimmax)])
            ylim([0 45])
            xlim([0 max(xlimmax)])
        end
    end
    xx = [xx; xlimmax];
    yy = [yy; ylimmax];

    F = getframe(f);
    imwrite(F.cdata, [cond_list{icond} '_distd_indiv_line_norm.png'], 'png')
    close
end

% avg
imp = [];
% cmap = [0.75 0 0; 0 0 0];
cmap = [0 0.75 0.75; 0 0 0];
xx = [];
yy = [];
for icond = 1:8
    f = figure(2100+icond); clf
    f.Position = [1176 85 602 846];
    ylimmax = [];
    xlimmax = [];

    for ilaser = 1:2

        subplot(2,1,ilaser)
        imp{icond,ilaser} = [];

        tmp = [];
        for isuccess = 1:2

            tmp_p = stack_distd{icond,ilaser}{isuccess};
            tmp = [];
            for itmp = 1:size(tmp_p,1)
                smp = tmp_p(itmp,:);
                
                if sum(isnan(smp)) == 1800
                    % do nothing
                else
                    fnan = isnan(smp);
                    if sum(fnan) < 1800-5
                        df_fnan = diff(fnan);
                        fdf_fnan = find(df_fnan == 1);
                        if length(fdf_fnan) == 1
                            final_pt = fdf_fnan;
                        else
                            final_pt = fdf_fnan(end);
                        end
                        smp = tmp_p(itmp,1:final_pt);
                        len_smp = length(smp);
                        x_smp = 1:((len_smp-1)/(1800-1)):(len_smp);
                        % x_smp = 1:(len_smp/(1800-1)):(len_smp+1);
                        y_interp = interp1(1:len_smp,smp,x_smp);
        
                        tmp = [tmp; y_interp];
                    end
                end
                
            end

            % tmp_avg = mean(stack_distd{icond,ilaser}{isuccess},'omitnan');
            tmp_avg = mean(tmp,'omitnan');
            % tmp_std = [];
            % for icol = 1:size(stack_distd{icond,ilaser}{isuccess},2)
            %     tmp_col = stack_distd{icond,ilaser}{isuccess}(:,icol);
            %     n_nan = isnan(tmp_col);
            %     n_nan = sum(n_nan);
            %     n_num = size(stack_distd{icond,ilaser}{isuccess},1)-n_nan;
            %     tmp_std(icol) = std(tmp_col,'omitnan')./sqrt(n_num);
            % end
            tmp_std = std(tmp,'omitnan')./sqrt(size(tmp,1));

            x_plot = 1:1800;
            y_plot = tmp_avg;
            e_plot = tmp_std;

            s_tmp = shadedErrorBar(x_plot,y_plot,e_plot); hold on
            s_tmp.mainLine.Color = cmap(isuccess,:);
            s_tmp.edge(1).LineStyle = 'none';
            s_tmp.edge(2).LineStyle = 'none';
            s_tmp.patch.FaceColor = cmap(isuccess,:);

            ax1 = gca;
            ylimmax = [ylimmax; ax1.YLim(2)];
            xlimmax = [xlimmax; ax1.XLim(2)];

            if ilaser == 1
                lasertext = 'Loff';
            elseif ilaser == 2
                lasertext = 'Lon';
            end

            % if isuccess == 2
            %     legend('fail','success')
            % end

            t_text = [cond_list{icond} ', ' lasertext];
            title(t_text)
        end
    end
    for ilaser = 1:2

        subplot(2,1,ilaser)
        for isuccess = 1:2

            % ylim([0 max(ylimmax)])
            ylim([0 45])
            xlim([0 max(xlimmax)])
        end
    end
    xx = [xx; xlimmax];
    yy = [ yy; ylimmax];

    F = getframe(f);
    imwrite(F.cdata, [cond_list{icond} '_distd_avg_line_norm.png'], 'png')
    close
end
















% cspd
imp = [];
for icond = 1:8
    f = figure(5000+icond); clf
    f.Position = [658 148 1205 745];
    ylimmax = [];
    xlimmax = [];
    climmax = [];

    stack_p = [];
    for ilaser = 1:2
        imp{icond,ilaser} = [];
        for isuccess = 1:2
            subplot(2,3,(ilaser-1)*3+isuccess), cla

            stack_p = [stack_p; stack_cspd{icond,ilaser}{isuccess}];

            tmp = mean(stack_cspd{icond,ilaser}{isuccess},'omitnan');
            fnan = isnan(tmp);
            f1 = find(fnan == 1);
            f1 = f1(1);

            tmp = stack_cspd{icond,ilaser}{isuccess}(:,1:f1);

            imp{icond,ilaser}{isuccess} = imagesc(tmp);

            ax1 = gca;
            ylimmax = [ylimmax; ax1.YLim(2)];
            xlimmax = [xlimmax; ax1.XLim(2)];
            climmax = [climmax; ax1.CLim(2)];

            if ilaser == 1
                lasertext = 'Loff';
            elseif ilaser == 2
                lasertext = 'Lon';
            end
            if isuccess == 1
                successtext = 'approach';
            elseif isuccess == 2
                successtext = 'intercept';
            end

            t_text = [cond_list{icond} ', ' successtext ', ' lasertext];
            title(t_text)

            colorbar
        end
    end
    for ilaser = 1:2
        for isuccess = 1:2
            subplot(2,3,(ilaser-1)*3+isuccess)

            ylim([0 max(ylimmax)])
            clim([0 max(climmax)])
        end
    end
    % stacked
    subplot(2,3,[3 6])
    imagesc(stack_p)
    xlim([0 max(xlimmax)])
    colorbar

    F = getframe(f);
    imwrite(F.cdata, [cond_list{icond} '_cspd_stacked_cmap.png'], 'png')
    close
end

% line chart:

imp = [];
% cmap = [0.75 0 0; 0 0 0];
cmap = [0 0.75 0.75; 0 0 0];
xx= [];
yy = [];
for icond = 1:8
    f = figure(3100+icond); clf
    f.Position = [1176 85 602 846];
    ylimmax = [];
    xlimmax = [];

    for ilaser = 1:2

        subplot(2,1,ilaser)
        imp{icond,ilaser} = [];
        for isuccess = 1:2

            tmp_p = stack_cspd{icond,ilaser}{isuccess};
            tmp = [];
            for itmp = 1:size(tmp_p,1)
                smp = tmp_p(itmp,:);
                
                if sum(isnan(smp)) == 1800
                    % do nothing
                else
                    fnan = isnan(smp);
                    if sum(fnan) < 1800-5
                        df_fnan = diff(fnan);
                        fdf_fnan = find(df_fnan == 1);
                        if length(fdf_fnan) == 1
                            final_pt = fdf_fnan;
                        else
                            final_pt = fdf_fnan(end);
                        end
                        smp = tmp_p(itmp,1:final_pt);
                        len_smp = length(smp);
                        x_smp = 1:((len_smp-1)/(1800-1)):(len_smp);
                        % x_smp = 1:(len_smp/(1800-1)):(len_smp+1);
                        y_interp = interp1(1:len_smp,smp,x_smp);
        
                        tmp = [tmp; y_interp];
                    end
                end
                
            end

            % tmp = mean(stack_cspd{icond,ilaser}{isuccess},'omitnan');
            tmp2 = plot(tmp'); hold on
            for itmp = 1:length(tmp2)
                % tmp(itmp).Color = [cmap(isuccess,:) 0.15];
                tmp2(itmp).Color = [cmap(isuccess,:) 0.2];
            end

            ax1 = gca;
            ylimmax = [ylimmax; ax1.YLim(2)];
            xlimmax = [xlimmax; ax1.XLim(2)];

            if ilaser == 1
                lasertext = 'Loff';
            elseif ilaser == 2
                lasertext = 'Lon';
            end

            % if isuccess == 2
            %     legend('fail','success')
            % end

            t_text = [cond_list{icond} ', ' lasertext];
            title(t_text)
        end
    end
    for ilaser = 1:2

        subplot(2,1,ilaser)
        for isuccess = 1:2

            % ylim([0 max(ylimmax)])
            ylim([0 50])
            xlim([0 max(xlimmax)])
        end
    end
    xx = [xx; xlimmax];
    yy = [yy ; ylimmax];

    % F = getframe(f);
    % imwrite(F.cdata, [cond_list{icond} '_cspd_indiv_line_norm.png'], 'png')
    % close
end

% avg
imp = [];
% cmap = [0.75 0 0; 0 0 0];
cmap = [0 0.75 0.75; 0 0 0];
xx = [];
yy = [];
for icond = 1:8
    f = figure(2100+icond); clf
    f.Position = [1176 85 602 846];
    ylimmax = [];
    xlimmax = [];

    for ilaser = 1:2

        subplot(2,1,ilaser)
        imp{icond,ilaser} = [];

        tmp = [];
        for isuccess = 1:2

            tmp_p = stack_cspd{icond,ilaser}{isuccess};
            tmp = [];
            for itmp = 1:size(tmp_p,1)
                smp = tmp_p(itmp,:);
                
                if sum(isnan(smp)) == 1800
                    % do nothing
                else
                    fnan = isnan(smp);
                    if sum(fnan) < 1800-5
                        df_fnan = diff(fnan);
                        fdf_fnan = find(df_fnan == 1);
                        if length(fdf_fnan) == 1
                            final_pt = fdf_fnan;
                        else
                            final_pt = fdf_fnan(end);
                        end
                        smp = tmp_p(itmp,1:final_pt);
                        len_smp = length(smp);
                        x_smp = 1:((len_smp-1)/(1800-1)):(len_smp);
                        % x_smp = 1:(len_smp/(1800-1)):(len_smp+1);
                        y_interp = interp1(1:len_smp,smp,x_smp);
        
                        tmp = [tmp; y_interp];
                    end
                end
                
            end

            % tmp_avg = mean(stack_cspd{icond,ilaser}{isuccess},'omitnan');
            tmp_avg = mean(tmp,'omitnan');
            % tmp_std = [];
            % for icol = 1:size(stack_cspd{icond,ilaser}{isuccess},2)
            %     tmp_col = stack_cspd{icond,ilaser}{isuccess}(:,icol);
            %     n_nan = isnan(tmp_col);
            %     n_nan = sum(n_nan);
            %     n_num = size(stack_cspd{icond,ilaser}{isuccess},1)-n_nan;
            %     tmp_std(icol) = std(tmp_col,'omitnan')./sqrt(n_num);
            % end
            tmp_std = std(tmp,'omitnan')./sqrt(size(tmp,1));

            x_plot = 1:1800;
            y_plot = tmp_avg;
            e_plot = tmp_std;

            s_tmp = shadedErrorBar(x_plot,y_plot,e_plot); hold on
            s_tmp.mainLine.Color = cmap(isuccess,:);
            s_tmp.edge(1).LineStyle = 'none';
            s_tmp.edge(2).LineStyle = 'none';
            s_tmp.patch.FaceColor = cmap(isuccess,:);

            ax1 = gca;
            ylimmax = [ylimmax; ax1.YLim(2)];
            xlimmax = [xlimmax; ax1.XLim(2)];

            if ilaser == 1
                lasertext = 'Loff';
            elseif ilaser == 2
                lasertext = 'Lon';
            end

            % if isuccess == 2
            %     legend('fail','success')
            % end

            t_text = [cond_list{icond} ', ' lasertext];
            title(t_text)
        end
    end
    for ilaser = 1:2

        subplot(2,1,ilaser)
        for isuccess = 1:2

            % ylim([0 max(ylimmax)])
            ylim([0 30])
            xlim([0 max(xlimmax)])
        end
    end
    xx = [xx; xlimmax];
    yy = [yy; ylimmax];

    F = getframe(f);
    imwrite(F.cdata, [cond_list{icond} '_cspd_avg_line_norm.png'], 'png')
    close
end








%%% testing to see how many of az are over 45?!

stat_dev = [];
for anim = 1:size(app_stat,2)
    for iapp = 1:size(app_stat{anim},2)

        tmp = find(app_stat{anim}{iapp}.az > 45);
        if isempty(tmp) == 0
            stat_dev = [stat_dev; anim iapp];
        end
    end
end

%%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%%
%%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%%
%%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%%
%%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%%
%%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%%
%%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%% - %%%

% distanced controlled plots for spd, az, and cspd
spd_prime = [];
az_prime = [];
cspd_prime = [];

spd_std = [];
az_std = [];
cspd_std = [];
for icond = 1:8
    for ilaser = 1:2
        spd_prime{icond,ilaser} = [];
        az_prime{icond,ilaser} = [];
        cspd_prime{icond,ilaser} = [];

        spd_std{icond,ilaser} = [];
        az_std{icond,ilaser} = [];
        cspd_std{icond,ilaser} = [];

        for isuccess = 1:2
            spd_prime{icond,ilaser}{isuccess} = [];
            az_prime{icond,ilaser}{isuccess} = [];
            cspd_prime{icond,ilaser}{isuccess} = [];

            spd_std{icond,ilaser}{isuccess} = [];
            az_std{icond,ilaser}{isuccess} = [];
            cspd_std{icond,ilaser}{isuccess} = [];
        end
    end
end

% x ranges
xcol = 0:5:80;

for icond = 1:8
    for ilaser = 1:2

        for isuccess = 1:2

            tmp_y1 = []; % spd
            tmp_y2 = []; % az
            tmp_y3 = []; % cspd
            for insert_col = 1:(length(xcol)-1)
                tmp_y1{insert_col} = []; % spd
                tmp_y2{insert_col} = []; % az
                tmp_y3{insert_col} = []; % cspd
            end


            for irow = 1:size(stack_distd{icond,ilaser}{isuccess},1)
                for icol = 1:size(stack_distd{icond,ilaser}{isuccess},2)

                    one_dat = stack_distd{icond,ilaser}{isuccess}(irow,icol);
                    if isnan(one_dat) == 0

                        xcol_min = find(xcol > one_dat);
                        if isempty(xcol_min) == 0
                            x_loc = min(xcol_min) - 1;
                        end
                        xcol_max = find(xcol < one_dat);
                        if isempty(xcol_max) == 0
                            x_loc = max(xcol_max);
                        end

                        tmp_spd = stack_spd{icond,ilaser}{isuccess}(irow,icol);
                        tmp_y1{x_loc} = [tmp_y1{x_loc}; tmp_spd];

                        tmp_az = stack_az{icond,ilaser}{isuccess}(irow,icol);
                        tmp_y2{x_loc} = [tmp_y2{x_loc}; tmp_az];

                        tmp_cspd = stack_cspd{icond,ilaser}{isuccess}(irow,icol);
                        tmp_y3{x_loc} = [tmp_y3{x_loc}; tmp_cspd];
                        
                    end
                end
            end

            for ip = 1:size(tmp_y1,2)

                spd_prime{icond,ilaser}{isuccess}(ip) = mean(tmp_y1{ip});
                az_prime{icond,ilaser}{isuccess}(ip) = mean(tmp_y2{ip});
                cspd_prime{icond,ilaser}{isuccess}(ip) = mean(tmp_y3{ip});

                spd_std{icond,ilaser}{isuccess}(ip) = std(tmp_y1{ip})/sqrt(size(tmp_y1{ip},1));
                az_std{icond,ilaser}{isuccess}(ip) = std(tmp_y2{ip})/sqrt(size(tmp_y2{ip},1));
                cspd_std{icond,ilaser}{isuccess}(ip) = std(tmp_y3{ip})/sqrt(size(tmp_y3{ip},1));
            end
        end
    end
end


x_plot = xcol-2.5;
x_plot(1) = [];
for icond = 1:8
    f = figure(1001+100*icond); clf
    f.Position = [718 498 1130 420];

    xlimmax = [];
    ylimmax = [];

    for ilaser = 1:2
        subplot(1,2,ilaser);
        for isuccess = 1:2

            s_tmp = shadedErrorBar(x_plot,spd_prime{icond,ilaser}{isuccess},spd_std{icond,ilaser}{isuccess}); hold on


            s_tmp.mainLine.Color = cmap(isuccess,:);
            s_tmp.edge(1).LineStyle = 'none';
            s_tmp.edge(2).LineStyle = 'none';
            s_tmp.patch.FaceColor = cmap(isuccess,:);

            ax1 = gca;
            xlimmax = [xlimmax; ax1.XLim(2)];
            ylimmax = [ylimmax; ax1.YLim(2)];
        end
    end

    for ilaser = 1:2
        subplot(1,2,ilaser)

        ax1 = gca;
        ax1.XLim(2) = max(xlimmax);
        ax1.YLim(2) = max(ylimmax);

        if ilaser == 1
            lasertext = 'Loff';
        elseif ilaser == 2
            lasertext = 'Lon';
        end
        t_text = [cond_list{icond} ', ' lasertext];
        title(t_text)

        ylabel('spd')
        xlabel('dist (cm)')
        
    end

    F = getframe(f);
    imwrite(F.cdata, [cond_list{icond} '_spd_distdctrl_cmap.png'], 'png')
    close
end

            
        
% az
x_plot = xcol-2.5;
x_plot(1) = [];
for icond = 1:8
    f = figure(2001+100*icond); clf
    f.Position = [718 498 1130 420];

    xlimmax = [];
    ylimmax = [];

    for ilaser = 1:2
        subplot(1,2,ilaser);
        for isuccess = 1:2

            s_tmp = shadedErrorBar(x_plot,az_prime{icond,ilaser}{isuccess},az_std{icond,ilaser}{isuccess}); hold on


            s_tmp.mainLine.Color = cmap(isuccess,:);
            s_tmp.edge(1).LineStyle = 'none';
            s_tmp.edge(2).LineStyle = 'none';
            s_tmp.patch.FaceColor = cmap(isuccess,:);

            ax1 = gca;
            xlimmax = [xlimmax; ax1.XLim(2)];
            ylimmax = [ylimmax; ax1.YLim(2)];
        end
    end

    for ilaser = 1:2
        subplot(1,2,ilaser)

        ax1 = gca;
        ax1.XLim(2) = max(xlimmax);
        ax1.YLim(2) = max(ylimmax);

        if ilaser == 1
            lasertext = 'Loff';
        elseif ilaser == 2
            lasertext = 'Lon';
        end
        t_text = [cond_list{icond} ', ' lasertext];
        title(t_text)

        ylabel('az')
        xlabel('dist (cm)')
        
    end

    F = getframe(f);
    imwrite(F.cdata, [cond_list{icond} '_az_distdctrl_cmap.png'], 'png')
    close
end




    
% cspd
x_plot = xcol-2.5;
x_plot(1) = [];
for icond = 1:8
    f = figure(2001+100*icond); clf
    f.Position = [718 498 1130 420];

    xlimmax = [];
    ylimmax = [];

    for ilaser = 1:2
        subplot(1,2,ilaser);
        for isuccess = 1:2

            s_tmp = shadedErrorBar(x_plot,cspd_prime{icond,ilaser}{isuccess},cspd_std{icond,ilaser}{isuccess}); hold on


            s_tmp.mainLine.Color = cmap(isuccess,:);
            s_tmp.edge(1).LineStyle = 'none';
            s_tmp.edge(2).LineStyle = 'none';
            s_tmp.patch.FaceColor = cmap(isuccess,:);

            ax1 = gca;
            xlimmax = [xlimmax; ax1.XLim(2)];
            ylimmax = [ylimmax; ax1.YLim(2)];
        end
    end

    for ilaser = 1:2
        subplot(1,2,ilaser)

        ax1 = gca;
        ax1.XLim(2) = max(xlimmax);
        ax1.YLim(2) = max(ylimmax);

        if ilaser == 1
            lasertext = 'Loff';
        elseif ilaser == 2
            lasertext = 'Lon';
        end
        t_text = [cond_list{icond} ', ' lasertext];
        title(t_text)

        ylabel('cspd')
        xlabel('dist (cm)')
        
    end

    F = getframe(f);
    imwrite(F.cdata, [cond_list{icond} '_cspd_distdctrl_cmap.png'], 'png')
    close
end