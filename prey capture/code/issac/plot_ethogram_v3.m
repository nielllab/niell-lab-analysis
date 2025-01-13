% capture rate during laser on vs laser off

clear
% close all

path1 = genpath('f:\prey_capture_lkhd_analysis\code\');
% addpath('f:\prey_capture_lkhd_analysis\code\')
addpath(path1);
% addpath('f:\code\')
path2 = genpath('f:\code\');
addpath(path2);
load('f:\prey_capture_lkhd_analysis\d_save\d_smp_5_dat+memo_v5.mat')
load('f:\prey_capture_lkhd_analysis\d_save\cdat_merged5_6_11_24.mat')
% load('f:\prey_capture_lkhd_analysis\d_save\app_d.mat')
% load('f:\prey_capture_lkhd_analysis\d_save\app_d_v4.mat')
% load('f:\prey_capture_lkhd_analysis\d_save\app_d_v5.mat')
% load('f:\prey_capture_lkhd_analysis\d_save\app_d_v11.mat')
load('f:\prey_capture_lkhd_analysis\d_save\app_d_v15.mat')
og_app_int = app_int;
load('f:\prey_capture_lkhd_analysis\d_save\exempt_list.mat')

% load 'capt_fr_out', 'out_dist', 'out_lkhd'
load('f:\prey_capture_lkhd_analysis\BinocOptoPreyCapture_out\pt995_thrs_extra.mat')

% load corrected time csv
f_name = 'csv_new_Def_v2.csv';
t_meta = readcell(['f:\prey_capture_lkhd_analysis\BinocOptoPreyCapture_out\' f_name],'DatetimeType','text');
t_meta_copy = t_meta;

csv_testing = readcell('f:\prey_capture_lkhd_analysis\csv_testing.csv','DatetimeType','text');

% rootdir = '/Volumes/T5 EVO/prey_capture_lkhd_analysis/';
rootdir = 'f:\prey_capture_lkhd_analysis\';

% load df_meta
f_name = 'df_meta.csv';
df_meta = readcell([rootdir f_name],'DatetimeType','text');

% cond_list = {'Wno','Wlb','Wsb','Lno','Lsb','Hsb','Wsw','Hno'};
cond_list = {'Wno','Lno','Hno','Wsw','Wlb','Wsb','Lsb','Hsb'};
% icond = 1;

original_df_meta = df_meta;

%%
% correct the t_meta (align it with capt_fr_out list)
for idf = 2:size(t_meta,1)
    disp(['Processing: ' num2mstr(idf) ' of ' num2str(size(t_meta,1))])
    if strcmp(t_meta{idf,1}(1),'P')
        f_anim = t_meta{idf,1};
        f_date = t_meta{idf,2};
        wheredash = find(f_date == '/');
        f_month = f_date(1:wheredash(1)-1);
        f_day = f_date((wheredash(1)+1):(wheredash(2)-1));
        if length(f_day) == 1
            f_day = ['0' f_day];
        end
        if length(f_month) == 1
            f_month = ['0' f_month];
        end
        f_date = [f_month f_day '21'];
        
        for itrial = 1:5
            atrial = itrial + 2;
    
            % if ismissing(t_meta{idf,atrial}) == 1
            if strcmp(t_meta{idf,atrial},'2nd last') == 0

                for ise = 1:size(capt_fr_out,1)
                    capt_name = capt_fr_out{ise,1};
                    whereund = find(capt_name == '_');
                    s_date = capt_name(1:whereund(1)-1);
                    s_anim = (capt_name(whereund(1)+1:whereund(2)-1));
                    s_trial = eval(capt_name(whereund(4)+1));
    
                    if strcmp(s_anim,f_anim) && strcmp(s_date,f_date) && s_trial == itrial
                        % disp(num2str(ise))
                        if isempty(capt_fr_out{ise,2}) == 0
                            t_meta{idf,atrial} = capt_fr_out{ise,2}(end);
                        else
                            t_meta{idf,atrial} = [];
                        end
                    end

                    if strcmp(csv_testing{idf,atrial+1}(1),'x')
                        t_meta{idf,atrial} = [];
                    elseif strcmp(csv_testing{idf,atrial+1}(1),'X')
                        t_meta{idf,atrial} = [];
                    end
                end

                if strcmp(t_meta_copy{idf,atrial}(1),'r')
                    % t_len = length(t_meta_copy{idf,atrial});
                    t_meta{idf,atrial} = eval(t_meta_copy{idf,atrial}(2:end));
                end
    
            elseif strcmp(t_meta{idf,atrial},'2nd last')
    
                for ise = 1:size(capt_fr_out,1)
                    capt_name = capt_fr_out{ise,1};
                    whereund = find(capt_name == '_');
                    s_date = capt_name(1:whereund(1)-1);
                    s_anim = (capt_name(whereund(1)+1:whereund(2)-1));
                    s_trial = eval(capt_name(whereund(4)+1));
    
                    if strcmp(s_anim,f_anim) && strcmp(s_date,f_date) && s_trial == itrial
                        t_meta{idf,atrial} = capt_fr_out{ise,2}(end-1);
                    end
                end

            end
    
        end
    end
end

%%
% correct captT first
for idf = 2:size(df_meta,1)

    id_name = df_meta{idf,17};
    id_date = num2str(df_meta{idf,18});
    if length(id_date) == 5
        id_month = id_date(1);
        id_day = id_date(2:3);
        id_year = id_date(4:5);
    elseif length(id_date) == 6
        id_month = id_date(1:2);
        id_day = id_date(3:4);
        if strcmp(id_day(1),'0')
            id_day = id_day(2);
        end
        id_year = id_date(5:6);
    end
    id_date = [id_month '/' id_day '/20' id_year];
    id_trial = df_meta{idf,16};
    id_cond = df_meta{idf,10};

    for f_itm = 2:size(t_meta,1)
        if strcmp(t_meta{f_itm,1},id_name) && strcmp(t_meta{f_itm,2},id_date)
            
            df_meta{idf,2} = t_meta{f_itm,id_trial+2}/60;

        end
    end

end

%% label each time frame as stationary, moving, approach, intercept, capture
csv_fldr = 'f:\prey_capture_lkhd_analysis\BinocOptoPreyCapture_wo_Spine\';

expctrltarget = 'Exp';
image_out = [];
cnt = 0;
match_id = [];

% cond_list = {'Wno','Wlb','Wsb','Lno','Lsb','Hsb','Wsw','Hno'};
targetcond = 8; %
% targetcond = 8; % Hsb
% targetcond = 6;
targetlaser = 'True';

app_d_loc = [];

for ifile = 2:size(df_meta,1)

    % identify which condition, laser, and ctrl/exp
    filecond = df_meta{ifile,9};
    expctrl = df_meta{ifile,10};
    laseryn = df_meta{ifile,13};

    % if strcmp(expctrl,expctrltarget) && strcmp(filecond,cond_list{targetcond}) 
    if strcmp(expctrl,expctrltarget) && strcmp(filecond,cond_list{targetcond}) && strcmp(laseryn,targetlaser)

        % ifile
        % laseryn = df_meta{ifile,13};
        if strcmp(laseryn,'True')
            lasert = 'y';
        elseif strcmp(laseryn,'False')
            lasert = 'n';
        end
        ftrial = df_meta{ifile,16};
        fname = df_meta{ifile,17};
        fdate = num2str(df_meta{ifile,18});
        if length(fdate) == 5
            fmonth = fdate(1);
            fday = fdate(2:3);
            if strcmp(fday(1),'0')
                fday = fday(2);
            end
            fyear = fdate(4:5);
        elseif length(fdate) == 6
            fmonth = fdate(1:2);
            fday = fdate(3:4);
            if strcmp(fday(1),'0')
                fday = fday(2);
            end
            fyear = fdate(5:6);
        end
        newfdate = [fmonth '/' fday '/20' fyear];
    
        condsav = [];
        for icond = 1:8
            if strcmp(cond_list{icond},filecond)
                condsav = icond;
            end
        end
        lasersav = [];
        if strcmp(laseryn,'False')
            lasersav = 1;
        elseif strcmp(laseryn,'True')
            lasersav = 2;
        end
    
        % print out for reading
        disp(['Processing: ' fdate ', ' fname ', ' num2str(ftrial) ' of 5 ------- ' num2str(ifile) ' of ' num2str(size(df_meta,1)) ', ' num2str(100*(ifile/size(df_meta,1))) '%'])
        
        capt_binary = [];
        capt_t = [];
        % use t_meta: (?)
        capt_dat = [];

        % use df_meta
        capt_dat = df_meta{ifile,2};
        capt_dat = round(capt_dat*60);
        
        if isempty(capt_dat) == 1
            aviname = [];
            if length(fdate) == 5
                fdate_p = ['0' fdate];
            else
                fdate_p = fdate;
            end
            if strcmp(laseryn,'True')
                flaser = 'y';
            else
                flaser = 'n';
            end
            fcond = filecond;
            ftrial_p = num2str(ftrial);
            aviname = [fdate_p '_' fname '_' flaser '_' fcond '_' ftrial_p '_TOP1.avi'];
        
            % load avi
            obj = VideoReader([csv_fldr aviname]);
            avi_frames = obj.NumFrames;
    
            capt_t = avi_frames;
            capt_binary = 0;
        else
            capt_t = capt_dat;
            if capt_t < 60*30
                capt_binary = 1;
            else
                capt_binary = 0; % still counts as failure if caught after 30s for this analysis
            end
        end


        % now work on all the rest of the parameters
        % find it in the app_h
        fapp = [];
        for napp = 1:size(app_h,2)
            if length(fdate) == 5
                app_fdate = ['0' fdate];
            else
                app_fdate = fdate;
            end
            if strcmp(app_fdate,app_h{napp}.fdate)
                if strcmp(fname,app_h{napp}.fname)
                    if strcmp(num2str(ftrial),app_h{napp}.ftrial)
                        fapp = napp;
                    end
                end
            end
        end
        % for v8, you need to add +1, for v9 you don't
        % fapp = fapp+1;
        fapp = fapp;

        app_d_loc = [app_d_loc; fapp];

        thrs_val = 60*30; % count only upto 1800 frames (30s mark)
        if capt_binary == 1 % this section trims the trial upto +60 fr from the time of capture
            if capt_t < thrs_val
                if capt_t + 60 < thrs_val
                    thrs_val = capt_t+60;
                end
            end
        end
        cnt_app = [];
        mark_app = [];
        if isempty(app_d{fapp}) == 0
            for iapp = 1:size(app_d{fapp},1)
                % mark as approach or intercept
                finddist = find(app_stat{fapp}{iapp}.distd<4);
                
                st = app_d{fapp}(iapp,1);
                ed = app_d{fapp}(iapp,2);
                ed1 = ed;
                if ed == size(app_p{fapp},2)
                    ed;
                else
                    ed = ed + 1;
                end
                if st < thrs_val && ed < thrs_val
                    cnt_app = [cnt_app; app_d{fapp}(iapp,:)];
                    if isempty(finddist) == 0
                        mark_app = [mark_app; 1];
                    else
                        mark_app = [mark_app; 0];
                    end
                elseif st < thrs_val && ed > thrs_val
                    % cnt_app = [cnt_app; app_d{fapp}(iapp,1) 60*30];
                    cnt_app = [cnt_app; app_d{fapp}(iapp,1) thrs_val];
                    if isempty(finddist) == 0
                        mark_app = [mark_app; 1];
                    else
                        mark_app = [mark_app; 0];
                    end
                end
            end
        end

        % annote cdat.d_dat{x}.param:
            % (1) frame_id
            % (2) dist_sk = speed of mouse
            % (3) az = azimuth
            % (4) dist_from_crkt (cm)
            % (5) dist_from_wall
            % (6) mid_head_x
            % (7) mid_head_y
            % (8) cricket_x
            % (9) cricket_y
            % (10) cricket_lkhd
            % (11) dist_crkt = speed of cricket

        % calculate time stationary
        % param_out = cdat.d_dat{cdat_id}.param;
        param_out = app_p{fapp};
        if size(param_out,2) > thrs_val
            spd_m = param_out(1,1:thrs_val);

            ang_m = param_out(2,1:thrs_val);
            % ang_m = abs(ang_m);

            dist_from_crkt = param_out(3,1:thrs_val);
        else
            spd_m = param_out(1,:);

            ang_m = param_out(2,:);

            dist_from_crkt = param_out(3,:);
        end

        % smooth it
        % spd_m = smooth(spd_m,30);
        % ang_m = smooth(ang_m,30);
        % dist_from_crkt = smooth(dist_from_crkt,30);
        % dist_from_wall = smooth(dist_from_wall,30);

        % create array where 0 == stationary(0-3cm/s), 1 = moving(3+ cm/s),
        % 2 == approach, 3== intercept, 4 == capture, NaN(?) == no trial
        len_fr = length(spd_m);

        tmp_nan_array = ones(1,1800);
        for ifr = 1:len_fr

            if spd_m(ifr) < 3
                tmp_nan_array(ifr) = 2; % stationary
            elseif spd_m(ifr) > 2.99999999999999
                tmp_nan_array(ifr) = 3; % mobile (3cm/s+)
            end

        end

        % approaches
        for iapp = 1:size(cnt_app,1)
            st = cnt_app(iapp,1);
            ed = cnt_app(iapp,2);

            for ifr = st:ed
                tmp_nan_array(ifr) = 4; % in approach
            end
        end

        % intercept (incoming, outgoing)
        for ifr = 2:len_fr
            if dist_from_crkt(ifr-1) > 4 % enter
                if dist_from_crkt(ifr) < 4
                    tmp_nan_array(ifr) = 5;
                end

            elseif dist_from_crkt(ifr-1) < 4 % exit
                if dist_from_crkt(ifr) > 4
                    tmp_nan_array(ifr) = 6;
                end
            end
        end

        % capture
        if capt_binary == 1
            tmp_nan_array(round(capt_t)) = 7;
        end

        cnt = cnt + 1;
        image_out = [image_out; tmp_nan_array];

        match_id = [match_id; ifile cnt fapp]; % useful for debugging, tracking down
    end
end


figg = figure(504); clf
figg.Position = [407 65 1499 900];
img1 = imagesc(image_out); hold on
colormap(jet);

% go through each line to put a visible marker
cmap = colormap(jet(7));
cmap(4,:) = [0 0 0]; colormap(cmap);
for iline = 1:size(image_out,1)

    tmp_array = image_out(iline,:);
    f5 = find(tmp_array == 5);
    if isempty(f5) == 0
        f6 = find(tmp_array == 6);
        for imarker = 1:length(f5)
            tmpp = plot(f5(imarker),iline,'.'); hold on
            tmpp.MarkerSize = 15;
            tmpp.Color = cmap(5,:);
            % connect with a line to the next one/leaving
            
            f6afterf5 = find(f6 > f5(imarker));
            if isempty(f6afterf5) == 0
                connectpt = f6(f6afterf5(1));
                tmpp = plot([f5(imarker) connectpt],[iline iline],'-','LineWidth',2); hold on
                tmpp.Color = cmap(5,:);
                f6(f6afterf5(1)) = [];
            end
        end
    end

    f6 = find(tmp_array == 6);
    if isempty(f6) == 0
        for imarker = 1:length(f6)
            tmpp = plot(f6(imarker),iline,'.'); hold on
            tmpp.MarkerSize = 15;
            tmpp.Color = cmap(6,:);
        end
    end

    f7 = find(tmp_array == 7);
    if isempty(f7) == 0
        tmpp = plot(f7,iline,'^'); hold on
        tmpp.MarkerSize = 8;
        tmpp.Color = cmap(7,:);
        tmpp.LineWidth = 2;
    end
    
end

colorbar

F = getframe(figg);
imwrite(F.cdata, ['ethnogram_' cond_list{targetcond} '_laser' targetlaser '.png'], 'png')
close

% re-arrange in order of capture time
new_image_out = [];

order_id = [];
for iarray = 1:size(image_out,1)
    tmp_array = image_out(iarray,:);
    f7 = find(tmp_array == 7);

    if isempty(f7) == 0
        order_id(iarray) = f7;
    else
        order_id(iarray) = 1800;
    end
end
[o_id,o_ev] = sort(order_id);

for iarray = 1:size(image_out,1)
    target_row = o_ev(iarray);
    new_image_out(iarray,:) = image_out(target_row,:);
end



figg = figure(518); clf
figg.Position = [407 65 1499 900];
img1 = imagesc(new_image_out); hold on
colormap(jet);

% go through each line to put a visible marker
cmap = colormap(jet(7));
cmap(4,:) = [0 0 0]; colormap(cmap);
for iline = 1:size(new_image_out,1)

    tmp_array = new_image_out(iline,:);
    f5 = find(tmp_array == 5);
    if isempty(f5) == 0
        f6 = find(tmp_array == 6);
        for imarker = 1:length(f5)
            tmpp = plot(f5(imarker),iline,'.'); hold on
            tmpp.MarkerSize = 15;
            tmpp.Color = cmap(5,:);
            % connect with a line to the next one/leaving
            
            f6afterf5 = find(f6 > f5(imarker));
            if isempty(f6afterf5) == 0
                connectpt = f6(f6afterf5(1));
                tmpp = plot([f5(imarker) connectpt],[iline iline],'-','LineWidth',2); hold on
                tmpp.Color = cmap(5,:);
                f6(f6afterf5(1)) = [];
            end
        end
    end

    f6 = find(tmp_array == 6);
    if isempty(f6) == 0
        for imarker = 1:length(f6)
            tmpp = plot(f6(imarker),iline,'.'); hold on
            tmpp.MarkerSize = 10;
            tmpp.Color = cmap(6,:);
        end
    end

    f7 = find(tmp_array == 7);
    if isempty(f7) == 0
        tmpp = plot(f7,iline,'^'); hold on
        tmpp.MarkerSize = 8;
        tmpp.Color = cmap(7,:);
        tmpp.LineWidth = 2;
    end
    
end

colorbar


F = getframe(figg);
imwrite(F.cdata, ['ethnogram_sorted_' cond_list{targetcond} '_laser' targetlaser '.png'], 'png')
close




%%



% look at the last time it enters the 4cm zone till capture
% keep the alignment (fastest catches first)

% catch = 7
% enter = 5
% exit = 6
new_image_out_2 = []; xlimcnt = [];
for iarray = 1:size(new_image_out,1)
    t_array = new_image_out(iarray,:);
    t_len = length(t_array);

    t1 = find(t_array == 1);
    t1_begin = [];
    if isempty(t1) == 0
        t1_begin = t1(1);
    end

    t7 = find(t_array == 7);
    if isempty(t7) == 0

        t5 = find(t_array == 5);
        t6 = find(t_array == 6); % not sure if needed

        t5_till = find(t5 < t7);
        t5_before_t7 = t5(t5_till(end));
        if isempty(t1) == 0
            if t5_before_t7 > 4 && t1_begin > t7+3
                insert_array = t_array((t5_before_t7-3):(t7+3));
    
                insert_len = length(insert_array);
                connect_ones = ones(1,(30*60)-insert_len);
    
                new_image_out_2(iarray,:) = [insert_array connect_ones];
                xlimcnt = [xlimcnt; insert_len iarray];
            else
                insert_array = t_array(t5_before_t7:t7);
                insert_array = [1 1 1 insert_array 1 1 1];
    
                insert_len = length(insert_array);
                connect_ones = ones(1,30*60-insert_len);
    
                new_image_out_2(iarray,:) = [insert_array connect_ones];
                xlimcnt = [xlimcnt; insert_len iarray];
            end
        elseif isempty(t1) == 1
            insert_array = t_array(t5_before_t7:t7);
            insert_array = [1 1 1 insert_array 1 1 1];

            insert_len = length(insert_array);
            connect_ones = ones(1,30*60-insert_len);

            new_image_out_2(iarray,:) = [insert_array connect_ones];
            xlimcnt = [xlimcnt; insert_len iarray];
        end
    elseif isempty(t7) == 1
        new_image_out_2(iarray,:) = new_image_out(iarray,:);
    end
end

figg = figure(538); clf
figg.Position = [407 65 1499 900];
img1 = imagesc(new_image_out_2); hold on
% go through each line to put a visible marker
cmap = colormap(jet(7));
cmap(4,:) = [0 0 0]; colormap(cmap);
for iline = 1:size(new_image_out_2,1)

    tmp_array = new_image_out_2(iline,:);
    f5 = find(tmp_array == 5);
    if isempty(f5) == 0
        for imarker = 1:length(f5)
            tmpp = plot(f5(imarker),iline,'.'); hold on
            tmpp.MarkerSize = 10;
            tmpp.Color = cmap(5,:);
        end
    end

    f6 = find(tmp_array == 6);
    if isempty(f6) == 0
        for imarker = 1:length(f6)
            tmpp = plot(f6(imarker),iline,'.'); hold on
            tmpp.MarkerSize = 10;
            tmpp.Color = cmap(6,:);
        end
    end

    f7 = find(tmp_array == 7);
    if isempty(f7) == 0
        tmpp = plot(f7,iline,'.'); hold on
        tmpp.MarkerSize = 10;
        tmpp.Color = cmap(7,:);
    end
    
end

% correct for xlim
xlim([-5 (max(xlimcnt(:,1))+5)])
colorbar