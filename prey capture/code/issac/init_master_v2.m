clear all
close all

addpath('f:\prey_capture_lkhd_analysis\code')  

% note: this code is a combination of compute_param and app_d ...
% first, it uses the csv to compute the parameters, saves the csv and the
% parameters, then uses the parameters to approaches and intercepts.
% biggest differences compared to prior code are (1) that it computes both
% parameters and approach related variables in parallel, so everything
% should be synced, and (2) smoothing and stitching/minmax are different.

% change:   (a) trial starts when cricket is properly registered by DLC
%           (give it -5 grace frames if available, as a means to buffer
%           jitter)
%           (b) trial ends +30s mark or when the cricket is caught,
%           whichever is faster
%           (c) don't analyze after cricket is caught (give it +5 grace
%           frames, as a means to buffer capture jitter)
%           (d) arena corners are derived from individual csv files, then
%           used to compute distance from walls
%           (e) 'distance from wall' will have 4 variables, one for each
%           wall: north, east, south, west

% -- load exp/ctrl csv -- use df_meta to analyze (this csv already weeds
% out which experiment to analyze and not analyze; source: Elliott)
cd('f:\prey_capture_lkhd_analysis')
f_name = 'df_meta.csv';
df_meta = readcell(f_name,'DatetimeType','text');

% -- load data
rootfldr = 'f:\prey_capture_lkhd_analysis\BinocOptoPreyCapture_wo_Spine';
cd(rootfldr)

% compute everything from csv's first within each loop
csv_list = dir('*.csv');
avi_list = dir('*.avi');

thrs = .99; % cutoff for lkhd
pxpercm = 14.5; % pixel per cm
framerate = 60;

wall_yes = 0; % using pre-defined wall information (if not, use original csv)
if wall_yes == 1
    % wall data
    ex_wall = readcell('f:\prey_capture_lkhd_analysis\d_save\ex_for_walls.csv','DatetimeType','text');
    % get left, right, top, bottom border coordinates
    ex_wall_label = ex_wall(2,:);
    tmp_ex_wall = ex_wall; tmp_ex_wall(1:3,:) = [];
    tmp_ex_wall = cell2mat(tmp_ex_wall);
    right = median(tmp_ex_wall(:,26)); %x-coord
    top = median(tmp_ex_wall(:,33)); %y-coord
    left = median(tmp_ex_wall(:,35)); %x-coord
    bottom = median(tmp_ex_wall(:,45)); %y-coord
end

% figure save fldr
% sav_fldr = 'f:\prey_capture_lkhd_analysis\BinocOptoPreyCapture_png_tmp\';
sav_fldr = 'F:\prey_capture_lkhd_analysis\BinocOptoPreyCapture_png_tmp\2024_12_13\';

% wall information fldr
wall_fldr = '\\Honey-Lake\Media\Data\BinocOptoPreyCapture_interp_wo_cricket';

%% calculate the time to capture, parameters, etc.

% param structure:
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

vid_out = 0;

% approach output
app_d = []; % outlines each approach
app_stat = []; % param stats for each approach epochs
app_h = []; % approach information (mouse id, date, trial number, condition, etc)
app_co = []; % approach conditions met or not (binary readouts of the 3 conditions)
app_p = []; % full parameter of trial
app_p_nan = []; % full parameter of trial, but leaves out confidence of cricket below 0.99 as NaN
app_p_wall = []; % parameters just for the wall information
app_int = []; % counts which approaches ended up as intercept (success)

% save output folders
png_fldr = 'f:\prey_capture_lkhd_analysis\d_save\tmp_png\';
avi_fldr = 'f:\prey_capture_lkhd_analysis\d_save\approach_avi_save_v11\';

prox_lmt = 4; % intercept condition

deg_lmt = 90/2; % angle condition
prc_lmt = 5; % speed condition
% innate negative distance (i.e. mouse closing in on the cricket)

% to confirm the lengths all match and lineup with the video files..
csv_match = [];
wcsv_match = [];

file_cnt = 0;
for ifile = 2:size(df_meta,1)
% for ifile = [49 50 51 52 382 383 384 385]
    file_cnt = file_cnt + 1;

    % load avi
    f_name = df_meta{ifile,17};
    f_date = df_meta{ifile,18};
    if length(num2str(f_date)) == 5
        f_date = ['0' num2str(f_date)];
    else
        f_date = num2str(f_date);
    end
    f_cond = df_meta{ifile,9};
    f_laser = df_meta{ifile,13};
    if strcmp(f_laser,'True')
        f_laser = 'y';
    elseif strcmp(f_laser,'False')
        f_laser = 'n';
    end
    f_trial = num2str(df_meta{ifile,16});
    f_full = [f_date '_' f_name '_' f_laser '_' f_cond '_' f_trial '_TOP1.avi'];

    app_h{file_cnt}.fname = f_name;
    app_h{file_cnt}.fdate = f_date;
    app_h{file_cnt}.fcond = f_cond;
    app_h{file_cnt}.ftrial = f_trial;
    app_h{file_cnt}.flaser = f_laser;

    % load video
    obj = VideoReader([rootfldr '\' f_full]);
    avi_frames = obj.NumFrames;
    vid = read(obj);

    % load csv
    c_full = [f_date '_' f_name '_' f_laser '_' f_cond '_' f_trial '_TOP1.csv'];
    csv1 = [];
    csv1 = readcell([rootfldr '\' c_full],'DatetimeType','text');
    
    if size(vid,4) == (size(csv1,1)-3)
        % do nothing
        % disp(['      --- csv and video lengths match --- '])
    else
        disp('      --- csv and video lengths DO NOT match!!! --- ')
        csv_match = [csv_match; ifile];
    end

    % load csv with wall info
    wcsv1 = [];
    w_full = [f_date '_' f_name '_' f_laser '_' f_cond '_' f_trial '_TOP1_interp.csv'];
    wcsv1 = readcell([wall_fldr '\' w_full],'DatetimeType','text');
    wcsv2 = convert1to2(wcsv1);

    if size(vid,4) == (size(wcsv1,1)-3)
        % do nothing
        % disp(['      --- csv and video lengths match --- '])
    else
        disp('      --- wcsv and video lengths DO NOT match!!! --- ')
        wcsv_match = [wcsv_match; ifile];
    end

    if isempty(wcsv1) == 0
        disp(['   Processing: ' f_full ', file #: ' num2str(file_cnt) ' of ' num2str(size(df_meta,1)-1) ', ' num2str(file_cnt/(size(df_meta,1)-1)*100) '%'])
    end

    % get average of the walls (use median function)
    wcsv2_med = median(wcsv2,1);
    % average out the two points per corner
    TR = [wcsv2_med(25)+wcsv2_med(28) wcsv2_med(26)+wcsv2_med(29)]/2;
    TL = [wcsv2_med(31)+wcsv2_med(34) wcsv2_med(32)+wcsv2_med(35)]/2;
    BL = [wcsv2_med(37)+wcsv2_med(40) wcsv2_med(38)+wcsv2_med(41)]/2;
    BR = [wcsv2_med(43)+wcsv2_med(46) wcsv2_med(44)+wcsv2_med(47)]/2;
    wcenter = (TR + TL + BL + BR)/4;
    wtop = (TR(2)+TL(2))/2; % y pos
    wbot = (BL(2)+BR(2))/2; % y pos
    wleft = (TL(1)+BL(1))/2; % x pos
    wright = (TR(1)+BR(1))/2; % x pos

    %% parameter calculations
    % csv2 = cell2mat(csv1(4:end,2:end));
    % updated:
    csv2 = convert1to2(csv1);

    % fill with NaNs
    for icol = 3:3:21
        for irow = 1:size(csv2,1)
            if csv2(irow,icol) < thrs
                csv2(irow,icol-2) = NaN;
                csv2(irow,icol-1) = NaN;
            end
        end
    end
    % smoothing
    for icol = 1:size(csv2,2)
        if length(find(icol == 3:3:21)) == 0
            inc = csv2(:,icol);
            onc = smooth(inc,15);
            csv2(:,icol) = onc;
        end
    end

    m_nose_xyz = csv2(:,1:3);
    m_lear_xyz = csv2(:,4:6);
    m_rear_xyz = csv2(:,7:9);
    m_bhead_xyz = csv2(:,10:12);
    m_tail_xyz = csv2(:,13:15);
    c1_xyz = csv2(:,16:18);
    c2_xyz = csv2(:,19:21);

    % get wall distances
    wdist_top = wtop-m_nose_xyz(:,2);
    wdist_bot = wbot-m_nose_xyz(:,2);
    wdist_left = wleft-m_nose_xyz(:,1);
    wdist_right = wright-m_nose_xyz(:,1);
    wdist_center = [wcenter(1)-m_nose_xyz(:,1) wcenter(2)-m_nose_xyz(:,2)];
    wdist = [];
    wdist.top = wdist_top;
    wdist.bot = wdist_bot;
    wdist.left = wdist_left;
    wdist.right = wdist_right;
    wdist.center = wdist_center;

    % (1) frame_id
    % (2) spd = speed of mouse
    % (3) az = azimuth
    % (4) dfc (cm)
    % (5) dfw (cm)
    % (6) mid_head_x
    % (7) mid_head_y
    % (8) cricket_x
    % (9) cricket_y
    % (10) cricket_lkhd
    % (11) dist_crkt = speed of cricket

    c_avg = [];
    for irow = 1:size(c1_xyz,1)
        if c1_xyz(irow,3) > thrs && c2_xyz(irow,3) > thrs
            c_avg(irow,:) = mean([c1_xyz(irow,:); c2_xyz(irow,:)]);
        elseif c1_xyz(irow,3) > thrs && c2_xyz(irow,3) < thrs
            c_avg(irow,:) = c1_xyz(irow,:);
        elseif c1_xyz(irow,3) < thrs && c2_xyz(irow,3) > thrs
            c_avg(irow,:) = c2_xyz(irow,:);
        elseif c1_xyz(irow,3) < thrs && c2_xyz(irow,3) < thrs
            c_avg(irow,:) = [NaN NaN mean([c1_xyz(irow,3) c2_xyz(irow,3)])];
        end
    end
    fnan = find(c_avg(:,3) < 0.99);

    % create 
    ear_avg_x = mean([m_lear_xyz(:,1) m_rear_xyz(:,1)],2);
    ear_avg_y = mean([m_lear_xyz(:,2) m_rear_xyz(:,2)],2);
    ear_avg = []; % clear any previously defined var
    ear_avg = [ear_avg_x ear_avg_y];

    % spd:
    % based on nose
    m_nose_xy = m_nose_xyz(:,1:2);
    m_nose_diff = diff(m_nose_xy);
    m_spd = sqrt((m_nose_diff(:,1).^2)+(m_nose_diff(:,2).^2));
    m_spd_15 = smooth(m_spd,15);
    spd = (m_spd_15*framerate)/pxpercm; % cm/s
    spd = [0; spd]; % add 0 since spd is a diff array
    % filter (for possible NaN)
    V = spd';
    X = ~isnan(V);
    Y = cumsum(X-diff([1,X])/2);
    Z = interp1(1:nnz(X),V(X),Y);
    param_spd = Z;
    param_spd_nan = param_spd;
    param_spd_nan(fnan) = NaN;

    % cricket spd:
    % based on "c_avg"
    c_avg_xy = c_avg(:,1:2);
    c_avg_diff = diff(c_avg_xy);
    c_spd = sqrt((c_avg_diff(:,1).^2)+(c_avg_diff(:,2).^2));
    c_spd_15 = smooth(c_spd,15);
    c_spd = (c_spd_15*framerate)/pxpercm; % cm/s
    c_spd = [0; c_spd]; % add 0 since spd is a diff array
    V = c_spd';
    X = ~isnan(V);
    Y = cumsum(X-diff([1,X])/2);
    Z = interp1(1:nnz(X),V(X),Y);
    param_c_spd = Z;
    param_c_spd_nan = param_c_spd;
    param_c_spd_nan(fnan) = NaN;

    % az:
    % azimuth based on the nose/center between ears and cricket
    % creat average of c1 and c2
    m_pt1 = m_nose_xyz(:,1:2); % m_traj = [m_pt1; m_pt2];
    m_pt2 = ear_avg(:,1:2);
    c_pt1 = c_avg(:,1:2); % c_traj = [c_pt1; c_pt2];
    % find angle from starting to ending pt
    x_diff = m_pt2-m_pt1; b = abs(x_diff);
    y_diff = m_pt2-m_pt1; a = abs(y_diff);
    hyp = sqrt(x_diff.^2 + y_diff.^2); c = hyp;   

    % nose to center-ear (a):
    side1 = m_nose_xyz(:,1)-ear_avg(:,1);
    side2 = m_nose_xyz(:,2)-ear_avg(:,2);
    a = sqrt(side1.^2 + side2.^2);
    % nose to cricket (b):
    side1 = c_avg(:,1)-m_nose_xyz(:,1);
    side2 = c_avg(:,2)-m_nose_xyz(:,2);
    b = sqrt(side1.^2 + side2.^2);
    % center-ear to cricket (c; hyp):
    side1 = ear_avg(:,1)-c_avg(:,1);
    side2 = ear_avg(:,2)-c_avg(:,2);
    c = sqrt(side1.^2 + side2.^2);
    % find angle C
    ang_c = [];
    for irow = 1:size(a,1)
        a1 = a(irow);
        b1 = b(irow);
        c1 = c(irow);
        if isnan(a1) == 0 && isnan(b1) == 0 && isnan(c1) == 0
            ang_c(irow,:) = acos((a1^2 + b1^2 - c1^2)/(2 * a1 * b1));
        else
            ang_c(irow,:) = NaN;
        end
    end
    % filter (for possible NaN)
    V = ang_c';
    X = ~isnan(V);
    Y = cumsum(X-diff([1,X])/2);
    Z = interp1(1:nnz(X),V(X),Y);
    param_ang = rad2deg(Z);
    param_ang = 180-param_ang;
    param_ang_nan = param_ang;
    param_ang_nan(fnan) = NaN;

    % disttocrkt:
    % use: c_avg, m_nose_xyz
    d2c = [];
    for iframe = 1:size(vid,4)
        if c_avg(iframe,3) > thrs && m_nose_xyz(iframe,3) > thrs
            xd = c_avg(iframe,1) - m_nose_xyz(iframe,1);
            yd = c_avg(iframe,2) - m_nose_xyz(iframe,2);
            xyd = sqrt(xd^2 + yd^2);
            d2c = [d2c xyd];
        else
            d2c = [d2c NaN];
        end
    end
    % filter (for possible NaN)
    V = d2c;
    X = ~isnan(V);
    Y = cumsum(X-diff([1,X])/2);
    Z = interp1(1:nnz(X),V(X),Y);
    param_d2c = Z/pxpercm;
    % param_d2c = smooth(param_d2c,15); % corrected from 30
    param_d2c = smooth(param_d2c,3);
    param_d2c_nan = param_d2c;
    param_d2c_nan(fnan) = NaN;

    % conditions for "approach" label:
    % 1. azimuth
    vid_len = size(vid,4);
    az_tmp = param_ang;
    % az_tmp = smooth(az_tmp,30);
    az_tmp = smooth(az_tmp,3);
    c2 = zeros(1,vid_len); % criterion 2
    f2 = find(az_tmp < deg_lmt); % 
    c2(f2) = 1;

    % 2. speed should be 5cm/s or above
    speed = param_spd;
    % speed = smooth(speed,30);
    c3 = zeros(1,vid_len);
    f3 = find(speed > prc_lmt);
    c3(f3) = 1;

    % 3. decrease in distance
    c5 = zeros(1,vid_len);
    d_mc = param_d2c;
    d_mc_diff = diff(d_mc);
    d_mc_diff = [0; d_mc_diff];
    d_mc_diff = smooth(d_mc_diff,30);
    f5 = find(d_mc_diff < 0);
    c5(f5) = 1;

    all_c = []; % erase any previous
    all_c = c2 + c3 + c5;

    % redefine all_c
    all_d = all_c;
    rep_c = find(all_c == 3);
    all_c = zeros(1,length(all_c));
    all_c(rep_c) = 1;

    cout_pckg = [c2; c3; c5; all_d];

    % count the number of approaches in the 'all_c' array
    approach_rec = [];
    cnt = 0;
    st_app = [];
    ed_app = [];
    for iframe = 1:(vid_len-1)
        
        % potential start of approach
        if all_c(iframe) == 0 && all_c(iframe+1) ~= 0
            st_app = iframe;
        end
        % delete it if st_app is 2nd last frame (it happens..)
        if st_app == vid_len-1
            st_app = [];
        end

        % if isempty(st_app) == 0
        if isempty(st_app+1) == 0 % it's the next one over that's the beginning of approach
            % if all_c(iframe) == 3
            % if all_c(iframe) == 1
            if all_c(iframe+1) == 1
                cnt = cnt + 1;
                % next_zero = find(all_c(iframe+1:end) == 0); % what is the '+1' for?
                % next_zero = find(all_c(iframe:end) == 0);
                next_zero = find(all_c(iframe+1:end) == 0);
                if isempty(next_zero) == 0
                    next_zero = next_zero(1);
                    ed_app = iframe + next_zero;
                    % approach_rec(cnt,:) = [st_app ed_app];
                    % approach_rec(cnt,:) = [st_app+1 ed_app];
                    approach_rec(cnt,:) = [st_app+1 ed_app-1];

                    st_app = [];
                    ed_app = [];
                elseif isempty(next_zero) == 1
                    next_zero = length(all_c);
                    ed_app = next_zero;
                    % approach_rec(cnt,:) = [st_app ed_app];
                    % approach_rec(cnt,:) = [st_app+1 ed_app];
                    approach_rec(cnt,:) = [st_app+1 ed_app-1];

                    st_app = [];
                    ed_app = [];
                end
            end
        end
    end
    % function for weeding out approaches that start <4cm, and also trim
    % approaches after first <4cm (intercept)
    p_spd = speed;
    p_az = az_tmp;
    p_d2c = d_mc;
    app_int_out = [];
    [approach_rec,app_int_out]=trim_approach(approach_rec,p_spd,p_az,p_d2c);

    og_approach_rec = approach_rec; % used for drawing

    % stitching:
    stitch_left_to_do = 1;
    n_approach_rec = approach_rec;
    n_app_int_out = app_int_out;
    % app_int_out_copy = app_int_out;
    while stitch_left_to_do == 1

        stitch_list = [];
        for iapp = 1:(size(n_approach_rec,1)-1)
            trial_diff = n_approach_rec(iapp+1,1) - n_approach_rec(iapp,2);
            % if trial_diff < 10 % 10 frames is about 170 ms
            if trial_diff < 15 % 15 frames is 250 ms
                stitch_list = [stitch_list; 1];
            else
                stitch_list = [stitch_list; 0];
            end
        end
        if sum(stitch_list) > 0
            tmp_approach_rec = [];
            new_app_int_out = [];
            for iapp = 1:(size(n_approach_rec,1)-1)
                if stitch_list(iapp,1) == 1
                    tmp_approach_rec = [tmp_approach_rec; n_approach_rec(iapp,1) n_approach_rec(iapp+1,2)];
                    %
                    new_app_int_out = [new_app_int_out; n_app_int_out(iapp) n_app_int_out(iapp+1)]; % this is the only time it gets stitched
                else
                    if iapp ~= 1
                        if stitch_list(iapp-1,1) == 0
                            tmp_approach_rec = [tmp_approach_rec; n_approach_rec(iapp,1) n_approach_rec(iapp,2)];
                            %
                            new_app_int_out = [new_app_int_out; n_app_int_out(iapp) n_app_int_out(iapp)]; % same info, just double for structure
                        end
                    elseif iapp == 1
                        tmp_approach_rec = [tmp_approach_rec; n_approach_rec(iapp,1:2)];
                        %
                        new_app_int_out = [new_app_int_out; n_app_int_out(iapp) n_app_int_out(iapp)]; % same info, just double for structure
                    end
                end
            end
            if stitch_list(end) == 0
                tmp_approach_rec = [tmp_approach_rec; n_approach_rec(end,1:2)];
                %
                new_app_int_out = [new_app_int_out; n_app_int_out(iapp) n_app_int_out(iapp)]; % same info, just double for structure

            end
            n_approach_rec = [];
            n_approach_rec = tmp_approach_rec;
            
            n_app_int_out = [];
            for iapp = 1:size(new_app_int_out,1)
                if sum(new_app_int_out(iapp,:)) > 0
                    n_app_int_out(iapp) = 1;
                else
                    n_app_int_out(iapp) = 0;
                end
            end

            % % fix the app_int_out
            % stitch_list_p = stitch_list;
            % stitch_yes = 0;
            % if sum(stitch_list_p) > 0
            %     stitch_yes = 1;
            % end
            % % 
            % while stitch_yes == 1
            % 
            %     n_app = length(stitch_list_p);
            %     if n_app > 1
            %         if stitch_list_p(1) == 0
            %             stitch_list_p(1) = [];
            %             new_app_int_out = [new_app_int_out app_int_out_copy(1)];
            %             app_int_out_copy(1) = [];
            %         else
            %             stitch_list_p(1) = [];
            % 
            %             if app_int_out_copy(1) == 0 && app_int_out_copy(2) == 0
            %                 new_app_int_out = [new_app_int_out 0];
            %             else
            %                 new_app_int_out = [new_app_int_out 1];
            %             end
            %             app_int_out_copy(1) = [];
            %         end
            %     else
            %         new_app_int_out = [new_app_int_out app_int_out_copy(1)];
            %         stitch_yes = 0;
            %     end
            % end

        else
            stitch_left_to_do = 0;
        end

    end
    approach_rec = [];
    approach_rec = n_approach_rec;
    %
    n_app_int_out;

    % take out all the ones less than 170 ms, about 10 frames
    tmp_approach = [];
    %
    tmp_app_int_out = [];
    for iapp = 1:size(approach_rec,1)
        st = approach_rec(iapp,1);
        ed = approach_rec(iapp,2);
        if ed-st > 10
            tmp_approach = [tmp_approach; approach_rec(iapp,:)];
            %
            tmp_app_int_out = [tmp_app_int_out; n_app_int_out(iapp)];
        end
    end
    approach_rec = [];
    approach_rec = tmp_approach;
    %
    app_int_out = tmp_app_int_out;

    % save the 3 parameters at the very least:
    stat_rec = [];
    for iapp = 1:size(approach_rec,1)
        stat_rec{iapp} = [];
    end
    for iapp = 1:size(approach_rec,1)
        st = approach_rec(iapp,1);
        ed = approach_rec(iapp,2);

        stat_rec{iapp}.speed = param_spd(st:ed);
        stat_rec{iapp}.az = param_ang(st:ed);
        stat_rec{iapp}.distd = param_d2c(st:ed);
        stat_rec{iapp}.cspd = param_c_spd(st:ed);
    end

    if vid_out == 1
        cd(sav_fldr)
        for iframe = 1:size(vid,4)
            f = figure(101); clf;
            f.Position = [237 178 1148 809];
            sub1 = subplot(4,2,[3 5]); cla
            imagesc(squeeze(vid(:,:,:,iframe))); hold on
            plot([m_nose_xyz(iframe,1) ear_avg(iframe,1)],[m_nose_xyz(iframe,2) ear_avg(iframe,2)],'y-','LineWidth',2); hold on % nose to ears
            plot([m_nose_xyz(iframe,1) c_avg(iframe,1)],[m_nose_xyz(iframe,2) c_avg(iframe,2)],'y-','LineWidth',2); hold on % nose to cricket
            plot([ear_avg(iframe,1) c_avg(iframe,1)],[ear_avg(iframe,2) c_avg(iframe,2)],'b--'); hold on
            ax1 = gca; ax1.YDir = 'normal';
            title({['frame # ' num2str(iframe) ' of ' num2str(size(vid,4))],[f_name ', ' f_date ', ' f_cond ', laser:' f_laser ', trial #' f_trial],['ang: ' num2str(param_ang(iframe))]})

            % approach labels
            if isempty(og_approach_rec) == 0
                
                sub2 = subplot(4,2,2); cla
                for iapp = 1:size(og_approach_rec,1)
                    plot([og_approach_rec(iapp,1) og_approach_rec(iapp,2)],[2 2],'b.-','LineWidth',1); hold on
                end
            end
            if isempty(approach_rec) == 0

                sub2 = subplot(4,2,2); 
                for iapp = 1:size(approach_rec,1)
                    plot([approach_rec(iapp,1) approach_rec(iapp,2)],[1 1],'b.-','LineWidth',1); hold on
                end

                xlim([0 size(vid,4)])
                ax1 = gca;
                ax1.XTick = [];
                ax1.YLabel.String = '1:sticheced app';
                ylim([0.75 2.25])
            end
            subplot(4,2,2);
            tmp_v = vline(iframe,'r'); hold on
            tmp_v.LineStyle = '--';
            

            % speed
            sub3 = subplot(4,2,4); cla
            % plot(smooth(param_spd,30),'k-','LineWidth',2); hold on
            plot(param_spd,'k-','LineWidth',4); hold on
            plot(param_spd_nan,'y-','LineWidth',2); hold on
            xlim([0 size(vid,4)])
            ax1 = gca;
            ax1.XTick = [];
            ax1.YLabel.String = 'speed (cm/s)';
            tmp_v = vline(iframe,'r'); hold on
            tmp_v.LineStyle = '--';
            title(['spd=' num2str(param_spd_nan(iframe))])


            % angle
            sub4 = subplot(4,2,6); cla
            % plot(smooth(param_ang,30),'k-','LineWidth',2); hold on
            plot(param_ang,'k-','LineWidth',4); hold on
            plot(param_ang_nan,'y-','LineWidth',2); hold on
            xlim([0 size(vid,4)])
            ylim([0 180])
            ax1 = gca;
            ax1.XTick = [];
            ax1.YLabel.String = 'angle (deg)';
            tmp_v = vline(iframe,'r'); hold on
            tmp_v.LineStyle = '--';
            title(['ang=' num2str(param_ang_nan(iframe))])


            % distance to cricket
            sub5 = subplot(4,2,8); cla
            % plot(param_d2c,'k-','LineWidth',2); hold on
            plot(param_d2c,'k-','LineWidth',4); hold on
            plot(param_d2c_nan,'y-','LineWidth',2); hold on
            xlim([0 size(vid,4)])
            ax1 = gca;
            ax1.YLabel.String = 'distance to crkt (cm)';
            tmp_v = vline(iframe,'r'); hold on
            tmp_v.LineStyle = '--';
            title(['d2c=' num2str(param_d2c_nan(iframe))])


            % frame number string generator
            if iframe < 10
                f_str = ['000' num2str(iframe)];
            elseif iframe > 9 && iframe < 100
                f_str = ['00' num2str(iframe)];
            elseif iframe > 99 && iframe < 1000
                f_str = ['0' num2str(iframe)];
            elseif iframe > 999 && iframe < 10000
                f_str = num2str(iframe);
            end

            F = getframe(f);
            imwrite(F.cdata, [ f_full(1:end-4) '_' f_str '.png'], 'png')

        end
        png_list = dir('*.png');
        writerObj = VideoWriter([f_full(1:end-4) '_angleTrace1_ammendd2c.avi']);
        writerObj.FrameRate = 60;
    
        % open the video writer
        open(writerObj);
        % write the frames to the video
        for ipng = 1:length(png_list)
            i_handle = imread([sav_fldr png_list(ipng).name]);
            % convert the image to a frame
            frame = im2frame(i_handle);
            writeVideo(writerObj, frame);
        end
        % close the writer object
        close(writerObj);

        delete([sav_fldr '\*.png'])
    end


    % save out and reset fldr
    cd(rootfldr)

    app_d{file_cnt} = approach_rec;
    app_stat{file_cnt} = stat_rec;
    app_co{file_cnt} = cout_pckg;
    app_p{file_cnt} = [param_spd; param_ang; param_d2c'; param_c_spd];
    app_p_nan{file_cnt} = [param_spd_nan; param_ang_nan; param_d2c_nan'; param_c_spd_nan];
    app_p_wall{file_cnt} = wdist; % wall distance parameters are here
    app_int{file_cnt} = app_int_out;

    
end


%% don't run:
% misc
if vid_out == 1
        cd(sav_fldr)
        for iframe = 1:size(vid,4)
            f = figure(101); clf;
            imagesc(squeeze(vid(:,:,:,iframe))); hold on
            plot([m_nose_xyz(iframe,1) ear_avg(iframe,1)],[m_nose_xyz(iframe,2) ear_avg(iframe,2)],'r-','LineWidth',2); hold on % nose to ears
            plot([m_nose_xyz(iframe,1) c_avg(iframe,1)],[m_nose_xyz(iframe,2) c_avg(iframe,2)],'r-','LineWidth',2); hold on % nose to cricket
            plot([ear_avg(iframe,1) c_avg(iframe,1)],[ear_avg(iframe,2) c_avg(iframe,2)],'b--'); hold on
            ax1 = gca; ax1.YDir = 'normal';
            title([num2str(rad2deg(ang_c(iframe)))])

            % frame number string generator
            if iframe < 10
                f_str = ['000' num2str(iframe)];
            elseif iframe > 9 && iframe < 100
                f_str = ['00' num2str(iframe)];
            elseif iframe > 99 && iframe < 1000
                f_str = ['0' num2str(iframe)];
            elseif iframe > 999 && iframe < 10000
                f_str = num2str(iframe);
            end

            F = getframe(f);
            imwrite(F.cdata, [ f_full(1:end-4) '_' f_str '.png'], 'png')

        end
        png_list = dir('*.png');
        writerObj = VideoWriter([f_full(1:end-4) '_angleTrace1.avi']);
        writerObj.FrameRate = 60;
    
        % open the video writer
        open(writerObj);
        % write the frames to the video
        for ipng = 1:length(png_list)
            i_handle = imread([sav_fldr png_list(ipng).name]);
            % convert the image to a frame
            frame = im2frame(i_handle);
            writeVideo(writerObj, frame);
        end
        % close the writer object
        close(writerObj);

        delete([sav_fldr '\*.png'])
    end