clear all
% same as 'plot_multiparam_v2_scatterv_successvsfail.m'
% but this is for approach vs intercepts (approach = fail, intercept =
% success)

%%

addpath('f:\prey_capture_lkhd_analysis\code')  
addpath('f:\code')  
addpath('f:\code\hline_vline\')
addpath('f:\code\povilaskarvelis\daboxplot\')
addpath('f:\code\densityScatterChart-1.2.0.0\')

% load 'app_d'
% load('f:\prey_capture_lkhd_analysis\d_save\app_d_v9.mat')
% load('f:\prey_capture_lkhd_analysis\d_save\app_d_v12.mat')
load('f:\prey_capture_lkhd_analysis\d_save\app_d_v15.mat')
% load corrected capture time 'df_meta'
load('f:\prey_capture_lkhd_analysis\d_save\df_meta_v1.mat')

% condition list
cond_list = {'Wno','Lno','Hno','Wsw','Wlb','Wsb','Lsb','Hsb'};

%
% this section compiles all the approach labeled frames and gathers the
% spd, az, and dcrkt for the given frame and puts it into 1 array for each
% condition and laser
tmspd_s = []; tmspd_f = [];
taz_s = []; taz_f = [];
tdcrkt_s = []; tdcrkt_f = [];
tcspd_s = []; tcspd_f = [];
for icond = 1:8
    for ilaser = 1:2
        tmspd_s{icond,ilaser} = []; tmspd_f{icond,ilaser} = [];
        taz_s{icond,ilaser} = []; taz_f{icond,ilaser} = [];
        tdcrkt_s{icond,ilaser} = []; tdcrkt_f{icond,ilaser} = [];
        tcspd_s{icond,ilaser} = []; tcspd_f{icond,ilaser} = [];
    end
end

% --- don't need trial avg? since dissecting it at approach level ---
% tmspd_trial_s = []; tmspd_trial_f = [];
% taz_trial_s = []; taz_trial_f = [];
% tdcrkt_trial_s = []; tdcrkt_trial_f = [];
% tcspd_s = []; tcspd_f = [];
% for icond = 1:8
%     for ilaser = 1:2
%         tmspd_trial_s{icond,ilaser} = []; tmspd_trial_f{icond,ilaser} = [];
%         taz_trial_s{icond,ilaser} = []; taz_trial_f{icond,ilaser} = [];
%         tdcrkt_trial_s{icond,ilaser} = []; tdcrkt_trial_f{icond,ilaser} = [];
%         tcspd_s{icond,ilaser} = []; tcspd_f{icond,ilaser} = [];
%     end
% end

capt_time2 = []; % another capt_time (previous one doesn't align for 1 of the conditions?!)
for icond = 1:8
    for ilaser = 1:2
        capt_time2{icond,ilaser} = [];
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
    % capt_binary = [];
    % capt_t = [];
    
    % find df_meta correspondence
    df_loc = [];
    capt_t = 'hello'; % it's a text so if it can't be saved because of the format, it'll error out
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

        file_id = ifile + 0; % correct for 'app_d_v8' here by +1

        disp(['Currently: ' num2str(ifile) ' of ' num2str(size(app_h,2)) ' ... ' fname ', ' fdate ', ' fcond ', ' flaser ', ' ftrial])

        % first, find the most recent 4cm zone entry from the capture time
        % then, drop the approaches that happen later than the above
        % approach; that one will be the one that is labeled 'successful
        % approach'

        capp_d = app_d{file_id};
        if isempty(capp_d) == 0
            if isempty(capt_t) == 0
    
                fex = [];
                fex = find(capp_d(:,1) > capt_t);
                if isempty(fex) == 0
                    capp_d(fex,:) = [];
                end
            end
        end

        % success (intercept) vs fail (no intercept)
        % 'sflist'
        if isempty(capp_d) == 0
          
            sflist = [];
            % for iapp = 1:size(capp_d,1)
            %     tmp_dcrkt = app_stat{file_id}{iapp}.distd;  
            %     tmpf = find(tmp_dcrkt < 4);
            %     if isempty(tmpf) == 0
            %         sflist(iapp,1) = 1;
            %     else
            %         sflist(iapp,1) = 0;
            %     end
            % end
            sflist = app_int{ifile};


            if isempty(capt_t) == 0
                cflist = []; % this indicates which of the approach led to the actual cathc
                for iapp = 1:size(capp_d,1)
                    st = capp_d(iapp,1);
                    ed = capp_d(iapp,2);
                    if st < capt_t+1 && ed > capt_t-1
                        cflist = [cflist; 1];
                    else
                        cflist = [cflist; 0];
                    end
    
                end
                % if you don't find it within the approach behaviors, you
                % have to assign it to the most recent/prior one
                if sum(cflist) == 0
                    tmpcapt = find(capp_d(:,2) < capt_t);
                    tmpcapt = tmpcapt(end);
                    cflist(tmpcapt) = 1;
                end
            else
                cflist = zeros(size(capp_d,1),1);
            end
        
        end


       




        % continue working (marked 3:58pm 11/25/2024)
        if isempty(capp_d) == 0 % isempty(app_d{file_id}) == 0

            tomean_spd = [];
            tomean_az = [];
            tomean_distd = [];

            for iapp = 1:size(capp_d,1) % size(app_d{file_id},1)
                % tmp_spd = app_stat{file_id}{iapp}.speed;
                % tmp_az = app_stat{file_id}{iapp}.az;
                % tmp_distd = app_stat{file_id}{iapp}.distd';

                st = capp_d(iapp,1);
                ed = capp_d(iapp,2);
                tmp_spd = app_p_nan{file_id}(1,st:ed);
                tmp_az = app_p_nan{file_id}(2,st:ed);
                tmp_distd = app_p_nan{file_id}(3,st:ed);
                tmp_cspd = app_p_nan{file_id}(4,st:ed);

                fzero = find(tmp_distd < 0);
                if isempty(fzero) == 0
                    tmp_distd(fzero) = [];
                    tmp_spd(fzero) = [];
                    tmp_az(fzero) = [];
                    tmp_cspd(fzero) = [];
                end

                mean_spd = mean(tmp_spd,'omitnan');
                mean_az = mean(tmp_az,'omitnan');
                mean_distd = mean(tmp_distd,'omitnan');
                mean_cspd = mean(tmp_cspd,'omitnan');

                sfd = sflist(iapp);
                if sfd == 0
                    tmspd_f{cond_loc,laser_loc} = [tmspd_f{cond_loc,laser_loc} mean_spd];
                    taz_f{cond_loc,laser_loc} = [taz_f{cond_loc,laser_loc} mean_az];
                    tdcrkt_f{cond_loc,laser_loc} = [tdcrkt_f{cond_loc,laser_loc} mean_distd];
                    tcspd_f{cond_loc,laser_loc} = [tcspd_f{cond_loc,laser_loc} mean_cspd];
                elseif sfd == 1
                    tmspd_s{cond_loc,laser_loc} = [tmspd_s{cond_loc,laser_loc} mean_spd];
                    taz_s{cond_loc,laser_loc} = [taz_s{cond_loc,laser_loc} mean_az];
                    tdcrkt_s{cond_loc,laser_loc} = [tdcrkt_s{cond_loc,laser_loc} mean_distd];
                    tcspd_s{cond_loc,laser_loc} = [tcspd_s{cond_loc,laser_loc} mean_cspd];
                end

            end

            % tmspd_trial{cond_loc,laser_loc} = [tmspd_trial{cond_loc,laser_loc}; mean(tomean_spd)];
            % taz_trial{cond_loc,laser_loc} = [taz_trial{cond_loc,laser_loc}; mean(tomean_az)];
            % tdcrkt_trial{cond_loc,laser_loc} = [tdcrkt_trial{cond_loc,laser_loc}; mean(tomean_distd)];




            if isempty(capt_t) == 1
                capt_time2{cond_loc,laser_loc} = [capt_time2{cond_loc,laser_loc}; NaN];
            else
                capt_time2{cond_loc,laser_loc} = [capt_time2{cond_loc,laser_loc}; capt_t];
            end
        end

    end
end

% tcspd tends to be fairly noisy due to poorer DLC estimations, so get rid
% of the NaN and negative values
tcspd_s_2 = [];
tcspd_f_2 = [];
for icond = 1:8
    for ilaser = 1:2
        f1 = find(tcspd_s{icond,ilaser} < 0);
        f2 = isnan(tcspd_s{icond,ilaser});
        f2 = find(f2 == 1);
        f3 = [f1 f2];
        f3 = unique(f3);

        tcspd_s_2{icond,ilaser} = tcspd_s{icond,ilaser};
        tcspd_s_2{icond,ilaser}(f3) = [];

        f1 = find(tcspd_f{icond,ilaser} < 0);
        f2 = isnan(tcspd_f{icond,ilaser});
        f2 = find(f2 == 1);
        f3 = [f1 f2];
        f3 = unique(f3);

        tcspd_f_2{icond,ilaser} = tcspd_f{icond,ilaser};
        tcspd_f_2{icond,ilaser}(f3) = [];
    end
end

cindex = [0 0.75 0.75; 0 0 0];

% cd('F:\prey_capture_lkhd_analysis\figure\2024_11_25')
cd('F:\prey_capture_lkhd_analysis\figure\2024_12_02')

% histogram plots
for iparam = 1:4 % 'tmspd_s', 'taz_s', 'tdcrkt_s', 'tcspd_s'

    if iparam == 1
        param_s = tmspd_s; param_f = tmspd_f;
    elseif iparam == 2
        param_s = taz_s; param_f = taz_f;
    elseif iparam == 3
        param_s = tdcrkt_s; param_f = tdcrkt_f;
    elseif iparam == 4
        param_s = tcspd_s_2; param_f = tcspd_f_2;
    end

    figh = figure(iparam*1000+icond*ilaser); clf
    % figh.Position = [1999 -65 776 1734];
    figh.Position = [1119 66 703 908];

    h1 = [];
    h2 = [];
    v1 = [];
    v2 = [];

    xlimcmp = [];
    ylimcmp = [];

    for icond = 1:8
        for ilaser = 1:2
            subplot(8,2,(icond-1)*2+ilaser); cla
            
            h1{icond,ilaser} = histogram(param_s{icond,ilaser}); hold on
            h2{icond,ilaser} = histogram(param_f{icond,ilaser}); hold on

            h1{icond,ilaser}.FaceColor = cindex(1,:);
            h2{icond,ilaser}.FaceColor = cindex(2,:);

            h1{icond,ilaser}.Normalization = 'probability';
            h2{icond,ilaser}.Normalization = 'probability';

            [ro,po] = ttest2(param_s{icond,ilaser},param_f{icond,ilaser});

            if po < 0.05
                if ilaser == 1
                    title(['*' cond_list{icond} ', Loff, po=' num2str(po)]);
                elseif ilaser == 2
                    title(['*' cond_list{icond} ', Lon, po=' num2str(po)]);
                end
            else
                if ilaser == 1
                    title([cond_list{icond} ', Loff, po=' num2str(po)]);
                elseif ilaser == 2
                    title([cond_list{icond} ', Lon, po=' num2str(po)]);
                end
            end

            ax1 = gca;
            xlimcmp = [xlimcmp; ax1.XLim];
            ylimcmp = [ylimcmp; ax1.YLim];

            
        end
    end

    % xinputmin = min(xlimcmp(:,1));
    % xinputmax = max(xlimcmp(:,2));
    xinputmin = 0;
    xinputmax = 50;

    xhistmin = floor(xinputmin);
    % xhistmax = ceil(xinputmax);
    xhistmin = 0;
    xhistmax = 50;

    yinputmin = min(ylimcmp(:,1));
    % yinputmax = max(ylimcmp(:,2));
    yinputmax = 0.3;


    for icond = 1:8
        for ilaser = 1:2

            subplot(8,2,(icond-1)*2+ilaser);

            ax1 = gca;
            ax1.XLim = [xinputmin xinputmax];
            ax1.YLim = [yinputmin yinputmax];

            h1{icond,ilaser}.BinLimits = [xhistmin xhistmax];
            h1{icond,ilaser}.NumBins = 20;

            h2{icond,ilaser}.BinLimits = [xhistmin xhistmax];
            h2{icond,ilaser}.NumBins = 20;

            h1_mean = mean(param_s{icond,ilaser},'omitnan');
            h2_mean = mean(param_f{icond,ilaser},'omitnan');

            v1{icond,ilaser} = vline(h1_mean); hold on
            v2{icond,ilaser} = vline(h2_mean); hold on

            v1{icond,ilaser}.Color = cindex(1,:);
            v1{icond,ilaser}.LineStyle = '-';
            v1{icond,ilaser}.LineWidth = 1;

            v2{icond,ilaser}.Color = cindex(2,:);
            v2{icond,ilaser}.LineStyle = '-';
            v2{icond,ilaser}.LineWidth = 1;

            if icond == 8
                if iparam == 1
                    xlabel('speed (cm/s)')
                elseif iparam == 2
                    xlabel('azimuth (deg)')
                elseif iparam == 3
                    xlabel('dist to crkt (cm)')
                elseif iparam == 4
                    xlabel('c speed (cm/s)')
                end
            end
        end
    end

    % F = getframe(figh);
    % if iparam == 1
    %     imwrite(F.cdata, 'spd.png', 'png')
    % elseif iparam == 2
    %     imwrite(F.cdata, 'az.png', 'png')
    % elseif iparam == 3
    %     imwrite(F.cdata, 'distcrkt.png', 'png')
    % elseif iparam == 4
    %     imwrite(F.cdata, 'cspd.png', 'png')
    % end
    % close


end

%% violin plots of above data
addpath('F:\code\povilaskarvelis\daviolinplot')
% cindex = [0 0.75 0.75; 0 0 0; 0 0.75 0.75; 0 0 0];
cindex = parula(4);
cindex(4,:) = [0 0 0];
c = cindex;
condition_names = cond_list;
group_inx = [1 2 3 4 5 6 7 8]; % not being used here
group_names = {'S Loff','S Lon','F Loff','F Lon'};

for iparam = 1:4

    if iparam == 1
        param_s = tmspd_s; param_f = tmspd_f;
    elseif iparam == 2
        param_s = taz_s; param_f = taz_f;
    elseif iparam == 3
        param_s = tdcrkt_s; param_f = tdcrkt_f;
    elseif iparam == 4
        param_s = tcspd_s_2; param_f = tcspd_f_2;
    end

    figj = figure(iparam*100+icond*ilaser); clf
    % figh.Position = [1999 -65 776 1734];
    figj.Position = [11 219 1886 540];

    % rearrange 'param_s' and 'param_f'
    s_laseroff_len = [];
    s_laseron_len = [];
    f_laseroff_len = [];
    f_laseron_len = [];
    for icond = 1:8
        s_laseroff_len = [s_laseroff_len; length(param_s{icond,1})];
        s_laseron_len = [s_laseron_len; length(param_s{icond,2})];

        f_laseroff_len = [f_laseroff_len; length(param_f{icond,1})];
        f_laseron_len = [f_laseron_len; length(param_f{icond,2})];
    end
    s_maxlaseroff_len = max(s_laseroff_len);
    s_maxlaseron_len = max(s_laseron_len);
    f_maxlaseroff_len = max(f_laseroff_len);
    f_maxlaseron_len = max(f_laseron_len);
    
    % laseroff s
    nan_p1 = nan(s_maxlaseroff_len,8);
    for icond = 1:8
        p1 = [];
        p1 = param_s{icond,1};
        p1 = p1';
        nan_p1(1:s_laseroff_len(icond),icond) = p1;
    end

    % laseron s
    nan_p2 = nan(s_maxlaseron_len,8);
    for icond = 1:8
        p2 = [];
        p2 = param_s{icond,2};
        p2 = p2';
        nan_p2(1:s_laseron_len(icond),icond) = p2;
    end

    % laseroff f
    nan_f1 = nan(f_maxlaseroff_len,8);
    for icond = 1:8
        f1 = [];
        f1 = param_f{icond,1};
        f1 = f1';
        nan_f1(1:f_laseroff_len(icond),icond) = f1;
    end

    % laseron f
    nan_f2 = nan(f_maxlaseron_len,8);
    for icond = 1:8
        f2 = [];
        f2 = param_f{icond,2};
        f2 = f2';
        nan_f2(1:f_laseron_len(icond),icond) = f2;
    end

    data3 = [];
    data3{1} = nan_p1;
    data3{2} = nan_p2;
    data3{3} = nan_f1;
    data3{4} = nan_f2;

    h = daviolinplot(data3,'groups',group_inx,'colors',c,'box',3,...
        'boxcolor','w','scatter',2,'jitter',1,'scattercolor','same',...
        'scattersize',10,'scatteralpha',0.7,'linkline',1,...
        'xtlabels', condition_names,...
        'legend',group_names);
    
    if iparam == 1
        title('mouse speed')
        ylabel('speed (cm/s)');
        ylim([-5 65]);
    elseif iparam == 2
        title('head angle/azimuth')
        ylabel('az (deg)')

    elseif iparam == 3
        title('distance between mouse and cricket')
        ylabel('dist (cm)')

    elseif iparam == 4
        title('cricket speed')
        ylabel('crkt spd (cm/s)')
        ylim([-5 50])
        
    end
    xlabel('conditions');

    xl = xlim; xlim([xl(1)-0.1, xl(2)+0.3]); % make more space for the legend
    set(h.sc,'MarkerEdgeColor','none');      % remove marker edge color
    set(gca,'FontSize',10);


    % F = getframe(figj);
    % if iparam == 1
    %     imwrite(F.cdata, 'v_spd.png', 'png')
    % elseif iparam == 2
    %     imwrite(F.cdata, 'v_az.png', 'png')
    % elseif iparam == 3
    %     imwrite(F.cdata, 'v_distcrkt.png', 'png')
    % elseif iparam == 4
    %     imwrite(F.cdata, 'v_cspd.png', 'png')
    % end
    % close
end





%%
% (a) looks at start of approach to 4cm
% (a1) classifies each approach as either fail or success
% (b) leaves out approaches that are already initially within 4cm (plots out
% the number/percentage of total approaches), or approaches that never
% reach 4cm (grabs the closest distance point) (how?)

% there are 2 versions of (a), one at absolute time (seconds) x-axis and
% another normalized from 0 to 1 (1 being 4cm enter zone)

targetexpctrl = 'Exp';
% targetexpctrl = 'Ctrl';

% version 'absolute'
app_s = []; 
app_f = [];
app_4cnt = []; % counts classification for each approaches
for icond = 1:8
    for ilaser = 1:2
        app_s{icond,ilaser} = [];
        app_f{icond,ilaser} = [];


        app_s{icond,ilaser}.spd = [];
        app_f{icond,ilaser}.spd = [];
        app_s{icond,ilaser}.az = [];
        app_f{icond,ilaser}.az = [];
        app_s{icond,ilaser}.dcrkt = [];
        app_f{icond,ilaser}.dcrkt = [];
        app_s{icond,ilaser}.cspd = [];
        app_f{icond,ilaser}.cspd = [];

        app_4cnt{icond,ilaser} = [];
    end
end
num_app = [];
for icond = 1:8
    for ilaser = 1:2
        num_app{icond,ilaser} = [];
    end
end
sanity_chk = [];

capt_dat = []; % binary (n/y; 0/1), time (max=1800)
for icond = 1:8
    for ilaser = 1:2
        capt_dat{icond,ilaser}  = [];
    end
end


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
    capt_t = 'hello'; % it's a text so if it can't be saved because of the format, it'll error out
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

    if isempty(capt_t) == 1
        capt_binary = 0;
    elseif isempty(capt_t) == 0
        if capt_t < 30*60 % it's in frames now
            capt_binary = 1;
        else
            capt_binary = 0;
        end
    end

    if strcmp(expctrl,targetexpctrl)
        if capt_binary == 0
            capt_dat{cond_loc,laser_loc} = [capt_dat{cond_loc,laser_loc}; capt_binary 1800];
        elseif capt_binary == 1
            capt_dat{cond_loc,laser_loc} = [capt_dat{cond_loc,laser_loc}; capt_binary capt_t];
        end
    end


    if strcmp(expctrl,targetexpctrl)

        file_id = ifile + 0; % correct for 'app_d_v8' here by +1

        disp(['Currently: ' num2str(ifile) ' of ' num2str(size(app_h,2)) ' ... ' fname ', ' fdate ', ' fcond ', ' flaser ', ' ftrial])

        % first, find the most recent 4cm zone entry from the capture time
        % then, drop the approaches that happen later than the above
        % approach; that one will be the one that is labeled 'successful
        % approach'

        capp_d = app_d{file_id};
        tmp_appd = capp_d;
        tmp_app_int = app_int{file_id};
        if isempty(tmp_appd) == 0
            % addition: 2025_01_10
            if capt_binary == 0 % no capture
                thrs = 30*60;
                exclude_list = find(tmp_appd(:,1) > thrs);
                if isempty(exclude_list) == 0
                    tmp_appd(exclude_list,:) = [];
                    tmp_app_int(exclude_list) = [];
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
                    tmp_app_int(exclude_list) = [];
                    if isempty(tmp_appd) == 0
                        if tmp_appd(end,2) > thrs
                            tmp_appd(end,2) = thrs;
                        end
                    end
                end
            end
            capp_d = [];
            capp_d = tmp_appd;
            % end addition
        end
        tmp_n_app = size(capp_d,1);
        num_app{cond_loc,laser_loc} = [num_app{cond_loc,laser_loc}; tmp_n_app];

        if isempty(capp_d) == 0
            
            n_app = size(capp_d,1);
            sp_distd = app_p_nan{file_id}(3,:);
            sp_distd = smooth(sp_distd,30); % need to smooth because some trials, there are NaN's around the 4cm zones
            len_sp_distd = length(sp_distd);
            c_line = []; % intercept/crossing the 4cm line or zone
            for ilen = 1:(len_sp_distd-1)
                if sp_distd(ilen) > 4 && sp_distd(ilen+1) < 4
                    c_line = [c_line; ilen];
                end
            end
            if isempty(c_line) == 1 % the smoothing function may have deteriorated the distd
                sp_distd = app_p_nan{file_id}(3,:);
                len_sp_distd = length(sp_distd);
                c_line = []; % intercept/crossing the 4cm line or zone
                for ilen = 1:(len_sp_distd-1)
                    if sp_distd(ilen) > 4 && sp_distd(ilen+1) < 4
                        c_line = [c_line; ilen];
                    end
                end
            end
        
            if isempty(capt_t) == 0
                [ar,br] = min(capt_t-c_line); % br is the index (ar is the diff value (in frames #))
                cr_fr = c_line(br);
    
                if n_app ~= 1
                    % conundrum:
                    % mouse may enter THEN approach, or approach THEN enter
                    % 4cm.. there is no correct behavior
                    % should i distinguish these 2 behaviors?

                    % for now:
                    % if there is only 1 approach, that is the approach
                    % that "leads" to the capture
                    % otherwise if the 4cm happens after the given
                    % approach, that approach is thrown out and the
                    % previous approach is used as the "final/success"
                    % approach
                    % ideally, the mouse does the approach behavior, then
                    % during the behavior or just prior to it, enters the
                    % 4cm zone, but in some cases, the mouse does enter the
                    % zone after the last approach behavior then captures
                    % it
                    % fcapp = find(capp_d(:,2) < cr_fr); % last frame of approach normally shouldn't be before the 4cm enter zone
                    fcapp = find(capp_d(:,1) > cr_fr);
                    if length(fcapp) == n_app
                        fcapp(end) = [];
                    end
                    capp_d(fcapp,:) = [];
                end

                sflist = zeros(1,size(capp_d,1));
                sflist(end) = 1;
                % overwrite:
                % sflist = app_int{ifile}';
                sflist = tmp_app_int;
            else
                sflist = zeros(1,size(capp_d,1));
                % overwrite:
                % sflist = app_int{ifile}';
                sflist = tmp_app_int;
            end
        end

        if isempty(capp_d) == 0
            % nan template:
            tmp_nan = nan(1,1000);

            % classify (1) started further than 4cm, reached 4cm, (2)
            % started further than 4cm, never reached 4cm, (3) started
            % within 4cm

            tmp_label = [];

            vid_len = size(app_p_nan{file_id},2);
            for iapp_prime = 1:size(capp_d,1)
                st = capp_d(iapp_prime,1);
                ed1 = capp_d(iapp_prime,2);
                ed = ed1;
                if ed1 == vid_len
                    ed;
                else
                    ed = ed + 1;
                end

                tmp_dcrkt = app_p_nan{file_id}(3,st:ed);
                f4 = find(tmp_dcrkt < 4); 
                p1 = tmp_dcrkt(1); % starting distance
                if isnan(p1)
                    p1 = tmp_dcrkt(2);
                end
                if p1 > 4 && isempty(f4) == 0
                    c_label = 1;
                    c_ins = [1 0 0];
                elseif p1 > 4 && isempty(f4) == 1
                    c_label = 2;
                    c_ins = [0 1 0];
                elseif p1 < 4
                    c_label = 3;
                    c_ins = [0 0 1];
                end
                
                tmp_label = [tmp_label; c_label c_ins];

                % app_4cnt{cond_loc,laser_loc} = [app_4cnt{cond_loc,laser_loc}; tmp_label];

                % if c_label == 1, grab the param time-series data
                % if c_label == 1

                % p_len = ed-st+1;
                % 
                % P_spd = app_p_nan{file_id}(1,st:ed);
                % P_az = app_p_nan{file_id}(2,st:ed);
                % P_dcrkt = app_p_nan{file_id}(3,st:ed);
                % P_cspd = app_p_nan{file_id}(4,st:ed);
                p_len = ed1-st+1;

                P_spd = app_p_nan{file_id}(1,st:ed1);
                P_az = app_p_nan{file_id}(2,st:ed1);
                P_dcrkt = app_p_nan{file_id}(3,st:ed1);
                P_cspd = app_p_nan{file_id}(4,st:ed1);

                % need to be same length (use the tmp_nan)
                t_spd = tmp_nan;
                t_spd(1:p_len) = P_spd;

                t_az = tmp_nan;
                t_az(1:p_len) = P_az;

                t_dcrkt = tmp_nan;
                t_dcrkt(1:p_len) = P_dcrkt;

                t_cspd = tmp_nan;
                t_cspd(1:p_len) = P_cspd;

                if sflist(iapp_prime) == 0 % fail
                    app_f{cond_loc,laser_loc}.spd = [app_f{cond_loc,laser_loc}.spd; t_spd];
                    app_f{cond_loc,laser_loc}.az = [app_f{cond_loc,laser_loc}.az; t_az];
                    app_f{cond_loc,laser_loc}.dcrkt = [app_f{cond_loc,laser_loc}.dcrkt; t_dcrkt];
                    app_f{cond_loc,laser_loc}.cspd = [app_f{cond_loc,laser_loc}.cspd; t_cspd];

                elseif sflist(iapp_prime) == 1 % success
                    app_s{cond_loc,laser_loc}.spd = [app_s{cond_loc,laser_loc}.spd; t_spd];
                    app_s{cond_loc,laser_loc}.az = [app_s{cond_loc,laser_loc}.az; t_az];
                    app_s{cond_loc,laser_loc}.dcrkt = [app_s{cond_loc,laser_loc}.dcrkt; t_dcrkt];
                    app_s{cond_loc,laser_loc}.cspd = [app_s{cond_loc,laser_loc}.cspd; t_cspd];

                end

                % end
            end

            app_4cnt{cond_loc,laser_loc} = [app_4cnt{cond_loc,laser_loc}; tmp_label];
            sanity_chk = [sanity_chk; size(capp_d,1) size(tmp_label,1) tmp_n_app];
        
        end

    end
end

% bar graph depicting #'s of different approach types (1, 2, 3)
f = figure(3844); clf
f.Position = [95 337 1759 673];

p_bar = [];
p_bar_prob = [];
for icond = 1:8
    for ilaser = 1:2

        f1 = find(app_4cnt{icond,ilaser}(:,1) == 1);
        f1 = length(f1);
        f2 = find(app_4cnt{icond,ilaser}(:,1) == 2);
        f2 = length(f2);
        f3 = find(app_4cnt{icond,ilaser}(:,1) == 3);
        f3 = length(f3);

        tcnt = f1+f2+f3;
        p1 = f1/tcnt;
        p2 = f2/tcnt;
        p3 = f3/tcnt;

        p_bar((icond-1)*2+ilaser,1:3) = [f1 f2 f3];
        p_bar_prob((icond-1)*2+ilaser,1:3) = [p1 p2 p3];
    end

end

subplot(2,1,1); cla
b1 = bar(p_bar);
ax1 = gca;
ax1.XTick = 1:16;
ax1.XTickLabel = {'Wno Loff','Wno Lon', ...
    'Lno Loff','Lno Lon', ...
    'Hno Loff','Hno Lon', ...
    'Wsw Loff','Wsw Lon', ...
    'Wlb Loff','Wlb Lon', ...
    'Wsb Loff','Wsb Lon', ...
    'Lsb Loff','Lsb Lon', ...
    'Hsb Loff','Hsb Lon'};
L1 = legend('start >4cm, reach <4cm','start >4cm, dont reach <4cm','start <4cm');
L1.Location = 'northwest';
ylabel('cnt (#)')
title('approach types and counts')

subplot(2,1,2); cla
b1 = bar(p_bar_prob);
ax1 = gca;
ax1.XTick = 1:16;
ax1.XTickLabel = {'Wno Loff','Wno Lon', ...
    'Lno Loff','Lno Lon', ...
    'Hno Loff','Hno Lon', ...
    'Wsw Loff','Wsw Lon', ...
    'Wlb Loff','Wlb Lon', ...
    'Wsb Loff','Wsb Lon', ...
    'Lsb Loff','Lsb Lon', ...
    'Hsb Loff','Hsb Lon'};
% L1 = legend('start >4cm, reach <4cm','start >4cm, dont reach <4cm','start <4cm');
% L1.Location = 'northwest';
ylim([0 1])
ylabel('prob (%/100)')
title('approach types and probability')

% without the 'starts <4cm'

% bar graph depicting #'s of different approach types (1, 2, 3)
f = figure(3845); clf
f.Position = [95 337 1759 673];

p_bar = [];
p_bar_prob = [];
for icond = 1:8
    for ilaser = 1:2

        f1 = find(app_4cnt{icond,ilaser}(:,1) == 1);
        f1 = length(f1);
        f2 = find(app_4cnt{icond,ilaser}(:,1) == 2);
        f2 = length(f2);
        % f3 = find(app_4cnt{icond,ilaser}(:,1) == 3);
        % f3 = length(f3);

        tcnt = f1+f2+f3;
        p1 = f1/tcnt;
        p2 = f2/tcnt;
        % p3 = f3/tcnt;

        % p_bar((icond-1)*2+ilaser,1:3) = [f1 f2 f3];
        p_bar((icond-1)*2+ilaser,1:2) = [f1 f2];
        % p_bar_prob((icond-1)*2+ilaser,1:3) = [p1 p2 p3];
        p_bar_prob((icond-1)*2+ilaser,1:2) = [p1 p2];
    end

end

subplot(2,1,1); cla
b1 = bar(p_bar);
ax1 = gca;
ax1.XTick = 1:16;
ax1.XTickLabel = {'Wno Loff','Wno Lon', ...
    'Lno Loff','Lno Lon', ...
    'Hno Loff','Hno Lon', ...
    'Wsw Loff','Wsw Lon', ...
    'Wlb Loff','Wlb Lon', ...
    'Wsb Loff','Wsb Lon', ...
    'Lsb Loff','Lsb Lon', ...
    'Hsb Loff','Hsb Lon'};
L1 = legend('start >4cm, reach <4cm','start >4cm, dont reach <4cm');
L1.Location = 'northwest';
ylabel('cnt (#)')
title('approach types and counts')

subplot(2,1,2); cla
b1 = bar(p_bar_prob);
ax1 = gca;
ax1.XTick = 1:16;
ax1.XTickLabel = {'Wno Loff','Wno Lon', ...
    'Lno Loff','Lno Lon', ...
    'Hno Loff','Hno Lon', ...
    'Wsw Loff','Wsw Lon', ...
    'Wlb Loff','Wlb Lon', ...
    'Wsb Loff','Wsb Lon', ...
    'Lsb Loff','Lsb Lon', ...
    'Hsb Loff','Hsb Lon'};
% L1 = legend('start >4cm, reach <4cm','start >4cm, dont reach <4cm','start <4cm');
% L1.Location = 'northwest';
ylim([0 1])
ylabel('prob (%/100)')
title('approach types and probability')

% figure; just for num_app
plt_num = [];
for icond = 1:8
    for ilaser = 1:2
        plt_num(icond,ilaser) = mean(num_app{icond,ilaser});
    end
end
plt_sum = [];
for icond = 1:8
    for ilaser = 1:2
        plt_sum(icond,ilaser) = sum(num_app{icond,ilaser});
    end
end


%% plot 'app_s' and 'app_f' which are the parameter data for approach->intercept only, segregated by success vs fail

addpath('F:\code\shadedErrorBar');

% spd
for icond = 1:8
    f = figure(4777+icond); clf
    f.Position = [53 97 1846 978];
    for ilaser = 1:2

        if ilaser == 1
            subplot(2,3,1); cla
            dat = app_f{icond,ilaser}.spd;
            plot(dat'); hold on
            plot(mean(dat,'omitnan'),'k-','LineWidth',2); hold on
            title('Loff, fail')
            dat1 = dat;
            xlim([0 200])

            subplot(2,3,2); cla
            dat = app_s{icond,ilaser}.spd;
            plot(dat'); hold on
            plot(mean(dat,'omitnan'),'k-','LineWidth',2); hold on
            title('Loff, success')
            dat2 = dat;
            xlim([0 200])

            subplot(2,3,3); cla
            tmp_y1 = mean(dat1,'omitnan');
            tmp_e1 = std(dat1,'omitnan');
            ecnt = [];
            for itick = 1:1000
                tmp_r = dat1(:,itick);
                tmp_isnan = isnan(tmp_r); tmp_isnan = 1-tmp_isnan;
                ecnt(itick) = sum(tmp_isnan);
            end
            tmp_e1e = tmp_e1./sqrt(ecnt);
            tmp_x1 = 1:1000;
            s1 = shadedErrorBar(tmp_x1,tmp_y1,tmp_e1e);
            hold on

            tmp_y1 = mean(dat2,'omitnan');
            tmp_e1 = std(dat2,'omitnan');
            ecnt = [];
            for itick = 1:1000
                tmp_r = dat2(:,itick);
                tmp_isnan = isnan(tmp_r); tmp_isnan = 1-tmp_isnan;
                ecnt(itick) = sum(tmp_isnan);
            end
            tmp_e1e = tmp_e1./sqrt(ecnt);
            tmp_x1 = 1:1000;
            s2 = shadedErrorBar(tmp_x1,tmp_y1,tmp_e1e);
            s2.patch.FaceColor = [0 0.75 0.75];
            xlim([0 200])

        elseif ilaser == 2

            subplot(2,3,4); cla
            dat = app_f{icond,ilaser}.spd;
            plot(dat'); hold on
            plot(mean(dat,'omitnan'),'k-','LineWidth',2); hold on
            title('Lon, fail')
            dat1 = dat;
            xlim([ 0 200])

            subplot(2,3,5); cla
            dat = app_s{icond,ilaser}.spd;
            plot(dat'); hold on
            plot(mean(dat,'omitnan'),'k-','LineWidth',2); hold on
            title('Lon, success')
            dat2 = dat;
            xlim([ 0 200])

            subplot(2,3,6); cla
            tmp_y1 = mean(dat1,'omitnan');
            tmp_e1 = std(dat1,'omitnan');
            ecnt = [];
            for itick = 1:1000
                tmp_r = dat1(:,itick);
                tmp_isnan = isnan(tmp_r); tmp_isnan = 1-tmp_isnan;
                ecnt(itick) = sum(tmp_isnan);
            end
            tmp_e1e = tmp_e1./sqrt(ecnt);
            tmp_x1 = 1:1000;
            s1 = shadedErrorBar(tmp_x1,tmp_y1,tmp_e1e);
            hold on

            tmp_y1 = mean(dat2,'omitnan');
            tmp_e1 = std(dat2,'omitnan');
            ecnt = [];
            for itick = 1:1000
                tmp_r = dat2(:,itick);
                tmp_isnan = isnan(tmp_r); tmp_isnan = 1-tmp_isnan;
                ecnt(itick) = sum(tmp_isnan);
            end
            tmp_e1e = tmp_e1./sqrt(ecnt);
            tmp_x1 = 1:1000;
            s2 = shadedErrorBar(tmp_x1,tmp_y1,tmp_e1e);
            s2.patch.FaceColor = [0 0.75 0.75];
            xlim([0 200])
        end
    end

    F = getframe(f);
    sav_title = [cond_list{icond} '_spd.png'];
    imwrite(F.cdata, sav_title, 'png')
    close
end

%
% az
for icond = 1:8
    f = figure(3777+icond); clf
    f.Position = [53 97 1846 978];
    for ilaser = 1:2

        if ilaser == 1
            subplot(2,3,1); cla
            dat = app_f{icond,ilaser}.az;
            plot(dat'); hold on
            plot(mean(dat,'omitnan'),'k-','LineWidth',2); hold on
            title('Loff, fail')
            dat1 = dat;
            xlim([0 200])

            subplot(2,3,2); cla
            dat = app_s{icond,ilaser}.az;
            plot(dat'); hold on
            plot(mean(dat,'omitnan'),'k-','LineWidth',2); hold on
            title('Loff, success')
            dat2 = dat;
            xlim([0 200])

            subplot(2,3,3); cla
            tmp_y1 = mean(dat1,'omitnan');
            tmp_e1 = std(dat1,'omitnan');
            ecnt = [];
            for itick = 1:1000
                tmp_r = dat1(:,itick);
                tmp_isnan = isnan(tmp_r); tmp_isnan = 1-tmp_isnan;
                ecnt(itick) = sum(tmp_isnan);
            end
            tmp_e1e = tmp_e1./sqrt(ecnt);
            tmp_x1 = 1:1000;
            s1 = shadedErrorBar(tmp_x1,tmp_y1,tmp_e1e);
            hold on

            tmp_y1 = mean(dat2,'omitnan');
            tmp_e1 = std(dat2,'omitnan');
            ecnt = [];
            for itick = 1:1000
                tmp_r = dat2(:,itick);
                tmp_isnan = isnan(tmp_r); tmp_isnan = 1-tmp_isnan;
                ecnt(itick) = sum(tmp_isnan);
            end
            tmp_e1e = tmp_e1./sqrt(ecnt);
            tmp_x1 = 1:1000;
            s2 = shadedErrorBar(tmp_x1,tmp_y1,tmp_e1e);
            s2.patch.FaceColor = [0 0.75 0.75];
            xlim([0 200])

        elseif ilaser == 2

            subplot(2,3,4); cla
            dat = app_f{icond,ilaser}.az;
            plot(dat'); hold on
            plot(mean(dat,'omitnan'),'k-','LineWidth',2); hold on
            title('Lon, fail')
            dat1 = dat;
            xlim([ 0 200])

            subplot(2,3,5); cla
            dat = app_s{icond,ilaser}.az;
            plot(dat'); hold on
            plot(mean(dat,'omitnan'),'k-','LineWidth',2); hold on
            title('Lon, success')
            dat2 = dat;
            xlim([ 0 200])

            subplot(2,3,6); cla
            tmp_y1 = mean(dat1,'omitnan');
            tmp_e1 = std(dat1,'omitnan');
            ecnt = [];
            for itick = 1:1000
                tmp_r = dat1(:,itick);
                tmp_isnan = isnan(tmp_r); tmp_isnan = 1-tmp_isnan;
                ecnt(itick) = sum(tmp_isnan);
            end
            tmp_e1e = tmp_e1./sqrt(ecnt);
            tmp_x1 = 1:1000;
            s1 = shadedErrorBar(tmp_x1,tmp_y1,tmp_e1e);
            hold on

            tmp_y1 = mean(dat2,'omitnan');
            tmp_e1 = std(dat2,'omitnan');
            ecnt = [];
            for itick = 1:1000
                tmp_r = dat2(:,itick);
                tmp_isnan = isnan(tmp_r); tmp_isnan = 1-tmp_isnan;
                ecnt(itick) = sum(tmp_isnan);
            end
            tmp_e1e = tmp_e1./sqrt(ecnt);
            tmp_x1 = 1:1000;
            s2 = shadedErrorBar(tmp_x1,tmp_y1,tmp_e1e);
            s2.patch.FaceColor = [0 0.75 0.75];
            xlim([0 200])
        end
    end

    F = getframe(f);
    sav_title = [cond_list{icond} '_az.png'];
    imwrite(F.cdata, sav_title, 'png')
    close
end



%
% dcrkt
for icond = 1:8
    f = figure(2777+icond); clf
    f.Position = [53 97 1846 978];
    for ilaser = 1:2

        if ilaser == 1
            subplot(2,3,1); cla
            dat = app_f{icond,ilaser}.dcrkt;
            plot(dat'); hold on
            plot(mean(dat,'omitnan'),'k-','LineWidth',2); hold on
            title('Loff, fail')
            dat1 = dat;
            xlim([0 200])

            subplot(2,3,2); cla
            dat = app_s{icond,ilaser}.dcrkt;
            plot(dat'); hold on
            plot(mean(dat,'omitnan'),'k-','LineWidth',2); hold on
            title('Loff, success')
            dat2 = dat;
            xlim([0 200])

            subplot(2,3,3); cla
            tmp_y1 = mean(dat1,'omitnan');
            tmp_e1 = std(dat1,'omitnan');
            ecnt = [];
            for itick = 1:1000
                tmp_r = dat1(:,itick);
                tmp_isnan = isnan(tmp_r); tmp_isnan = 1-tmp_isnan;
                ecnt(itick) = sum(tmp_isnan);
            end
            tmp_e1e = tmp_e1./sqrt(ecnt);
            tmp_x1 = 1:1000;
            s1 = shadedErrorBar(tmp_x1,tmp_y1,tmp_e1e);
            hold on

            tmp_y1 = mean(dat2,'omitnan');
            tmp_e1 = std(dat2,'omitnan');
            ecnt = [];
            for itick = 1:1000
                tmp_r = dat2(:,itick);
                tmp_isnan = isnan(tmp_r); tmp_isnan = 1-tmp_isnan;
                ecnt(itick) = sum(tmp_isnan);
            end
            tmp_e1e = tmp_e1./sqrt(ecnt);
            tmp_x1 = 1:1000;
            s2 = shadedErrorBar(tmp_x1,tmp_y1,tmp_e1e);
            s2.patch.FaceColor = [0 0.75 0.75];
            xlim([0 200])

        elseif ilaser == 2

            subplot(2,3,4); cla
            dat = app_f{icond,ilaser}.dcrkt;
            plot(dat'); hold on
            plot(mean(dat,'omitnan'),'k-','LineWidth',2); hold on
            title('Lon, fail')
            dat1 = dat;
            xlim([ 0 200])

            subplot(2,3,5); cla
            dat = app_s{icond,ilaser}.dcrkt;
            plot(dat'); hold on
            plot(mean(dat,'omitnan'),'k-','LineWidth',2); hold on
            title('Lon, success')
            dat2 = dat;
            xlim([ 0 200])

            subplot(2,3,6); cla
            tmp_y1 = mean(dat1,'omitnan');
            tmp_e1 = std(dat1,'omitnan');
            ecnt = [];
            for itick = 1:1000
                tmp_r = dat1(:,itick);
                tmp_isnan = isnan(tmp_r); tmp_isnan = 1-tmp_isnan;
                ecnt(itick) = sum(tmp_isnan);
            end
            tmp_e1e = tmp_e1./sqrt(ecnt);
            tmp_x1 = 1:1000;
            s1 = shadedErrorBar(tmp_x1,tmp_y1,tmp_e1e);
            hold on

            tmp_y1 = mean(dat2,'omitnan');
            tmp_e1 = std(dat2,'omitnan');
            ecnt = [];
            for itick = 1:1000
                tmp_r = dat2(:,itick);
                tmp_isnan = isnan(tmp_r); tmp_isnan = 1-tmp_isnan;
                ecnt(itick) = sum(tmp_isnan);
            end
            tmp_e1e = tmp_e1./sqrt(ecnt);
            tmp_x1 = 1:1000;
            s2 = shadedErrorBar(tmp_x1,tmp_y1,tmp_e1e);
            s2.patch.FaceColor = [0 0.75 0.75];
            xlim([0 200])
        end
    end

    F = getframe(f);
    sav_title = [cond_list{icond} '_dcrkt.png'];
    imwrite(F.cdata, sav_title, 'png')
    close
end


%
% cspd
for icond = 1:8
    f = figure(3777+icond); clf
    f.Position = [53 97 1846 978];
    for ilaser = 1:2

        if ilaser == 1
            subplot(2,3,1); cla
            dat = app_f{icond,ilaser}.cspd;
            plot(dat'); hold on
            plot(mean(dat,'omitnan'),'k-','LineWidth',2); hold on
            title('Loff, fail')
            dat1 = dat;
            xlim([0 200])

            subplot(2,3,2); cla
            dat = app_s{icond,ilaser}.cspd;
            plot(dat'); hold on
            plot(mean(dat,'omitnan'),'k-','LineWidth',2); hold on
            title('Loff, success')
            dat2 = dat;
            xlim([0 200])

            subplot(2,3,3); cla
            tmp_y1 = mean(dat1,'omitnan');
            tmp_e1 = std(dat1,'omitnan');
            ecnt = [];
            for itick = 1:1000
                tmp_r = dat1(:,itick);
                tmp_isnan = isnan(tmp_r); tmp_isnan = 1-tmp_isnan;
                ecnt(itick) = sum(tmp_isnan);
            end
            tmp_e1e = tmp_e1./sqrt(ecnt);
            tmp_x1 = 1:1000;
            s1 = shadedErrorBar(tmp_x1,tmp_y1,tmp_e1e);
            hold on

            tmp_y1 = mean(dat2,'omitnan');
            tmp_e1 = std(dat2,'omitnan');
            ecnt = [];
            for itick = 1:1000
                tmp_r = dat2(:,itick);
                tmp_isnan = isnan(tmp_r); tmp_isnan = 1-tmp_isnan;
                ecnt(itick) = sum(tmp_isnan);
            end
            tmp_e1e = tmp_e1./sqrt(ecnt);
            tmp_x1 = 1:1000;
            s2 = shadedErrorBar(tmp_x1,tmp_y1,tmp_e1e);
            s2.patch.FaceColor = [0 0.75 0.75];
            xlim([0 200])

        elseif ilaser == 2

            subplot(2,3,4); cla
            dat = app_f{icond,ilaser}.cspd;
            plot(dat'); hold on
            plot(mean(dat,'omitnan'),'k-','LineWidth',2); hold on
            title('Lon, fail')
            dat1 = dat;
            xlim([ 0 200])

            subplot(2,3,5); cla
            dat = app_s{icond,ilaser}.cspd;
            plot(dat'); hold on
            plot(mean(dat,'omitnan'),'k-','LineWidth',2); hold on
            title('Lon, success')
            dat2 = dat;
            xlim([ 0 200])

            subplot(2,3,6); cla
            tmp_y1 = mean(dat1,'omitnan');
            tmp_e1 = std(dat1,'omitnan');
            ecnt = [];
            for itick = 1:1000
                tmp_r = dat1(:,itick);
                tmp_isnan = isnan(tmp_r); tmp_isnan = 1-tmp_isnan;
                ecnt(itick) = sum(tmp_isnan);
            end
            tmp_e1e = tmp_e1./sqrt(ecnt);
            tmp_x1 = 1:1000;
            s1 = shadedErrorBar(tmp_x1,tmp_y1,tmp_e1e);
            hold on

            tmp_y1 = mean(dat2,'omitnan');
            tmp_e1 = std(dat2,'omitnan');
            ecnt = [];
            for itick = 1:1000
                tmp_r = dat2(:,itick);
                tmp_isnan = isnan(tmp_r); tmp_isnan = 1-tmp_isnan;
                ecnt(itick) = sum(tmp_isnan);
            end
            tmp_e1e = tmp_e1./sqrt(ecnt);
            tmp_x1 = 1:1000;
            s2 = shadedErrorBar(tmp_x1,tmp_y1,tmp_e1e);
            s2.patch.FaceColor = [0 0.75 0.75];
            xlim([0 200])
        end
    end

    F = getframe(f);
    sav_title = [cond_list{icond} '_cspd.png'];
    imwrite(F.cdata, sav_title, 'png')
    close
end



%%

% now as a time-series % (0 to 1)
interp_f = [];
interp_s = [];
for icond = 1:8
    for ilaser = 1:2
        interp_f{icond,ilaser} = [];
        interp_s{icond,ilaser} = [];


        interp_f{icond,ilaser}.spd = [];
        interp_s{icond,ilaser}.spd = [];
        interp_f{icond,ilaser}.az = [];
        interp_s{icond,ilaser}.az = [];
        interp_f{icond,ilaser}.dcrkt = [];
        interp_s{icond,ilaser}.dcrkt = [];
        interp_f{icond,ilaser}.cspd = [];
        interp_s{icond,ilaser}.cspd = [];
    end
end

smp = 100;
for icond =  1:8
    for ilaser = 1:2

        n_app = size(app_s{icond,ilaser}.spd,1);
        for iapp = 1:n_app
            % spd
            tmp = app_s{icond,ilaser}.spd(iapp,:);
            tmp_nan = isnan(tmp);
            tmp(tmp_nan) = [];
            tmp_len = length(tmp);

            if tmp_len > 5
                tmp_x = 1:tmp_len;
                tmp_y = tmp;
                tmp_interp = 1:((tmp_len-1)/smp):tmp_len;
    
                vq = interpn(tmp_x,tmp_y,tmp_interp,'spline');
    
                interp_s{icond,ilaser}.spd = [interp_s{icond,ilaser}.spd; vq];
            end
        end

        n_app = size(app_s{icond,ilaser}.az,1);
        for iapp = 1:n_app
            % az
            tmp = app_s{icond,ilaser}.az(iapp,:);
            tmp_nan = isnan(tmp);
            tmp(tmp_nan) = [];
            tmp_len = length(tmp);

            if tmp_len > 5
                tmp_x = 1:tmp_len;
                tmp_y = tmp;
                tmp_interp = 1:((tmp_len-1)/smp):tmp_len;
    
                vq = interpn(tmp_x,tmp_y,tmp_interp,'spline');
    
                interp_s{icond,ilaser}.az = [interp_s{icond,ilaser}.az; vq];
            end
        end

        n_app = size(app_s{icond,ilaser}.dcrkt,1);
        for iapp = 1:n_app
            % dcrkt
            tmp = app_s{icond,ilaser}.dcrkt(iapp,:);
            tmp_nan = isnan(tmp);
            tmp(tmp_nan) = [];
            tmp_len = length(tmp);

            if tmp_len > 5
                tmp_x = 1:tmp_len;
                tmp_y = tmp;
                tmp_interp = 1:((tmp_len-1)/smp):tmp_len;
    
                vq = interpn(tmp_x,tmp_y,tmp_interp,'spline');
    
                interp_s{icond,ilaser}.dcrkt = [interp_s{icond,ilaser}.dcrkt; vq];
            end
        end

        n_app = size(app_s{icond,ilaser}.cspd,1);
        for iapp = 1:n_app
            % cspd
            tmp = app_s{icond,ilaser}.cspd(iapp,:);
            tmp_nan = isnan(tmp);
            tmp(tmp_nan) = [];
            tmp_len = length(tmp);

            if tmp_len > 5
                tmp_x = 1:tmp_len;
                tmp_y = tmp;
                tmp_interp = 1:((tmp_len-1)/smp):tmp_len;
    
                vq = interpn(tmp_x,tmp_y,tmp_interp,'spline');
    
                interp_s{icond,ilaser}.cspd = [interp_s{icond,ilaser}.cspd; vq];
            end
        end



        % f



        n_app = size(app_f{icond,ilaser}.spd,1);
        for iapp = 1:n_app
            % spd
            tmp = app_f{icond,ilaser}.spd(iapp,:);
            tmp_nan = isnan(tmp);
            tmp(tmp_nan) = [];
            tmp_len = length(tmp);

            if tmp_len > 5
                tmp_x = 1:tmp_len;
                tmp_y = tmp;
                tmp_interp = 1:((tmp_len-1)/smp):tmp_len;
    
                vq = interpn(tmp_x,tmp_y,tmp_interp,'spline');
    
                interp_f{icond,ilaser}.spd = [interp_f{icond,ilaser}.spd; vq];
            end
        end

        n_app = size(app_f{icond,ilaser}.az,1);
        for iapp = 1:n_app
            % az
            tmp = app_f{icond,ilaser}.az(iapp,:);
            tmp_nan = isnan(tmp);
            tmp(tmp_nan) = [];
            tmp_len = length(tmp);

            if tmp_len > 5
                tmp_x = 1:tmp_len;
                tmp_y = tmp;
                tmp_interp = 1:((tmp_len-1)/smp):tmp_len;
    
                vq = interpn(tmp_x,tmp_y,tmp_interp,'spline');
    
                interp_f{icond,ilaser}.az = [interp_f{icond,ilaser}.az; vq];
            end
        end

        n_app = size(app_f{icond,ilaser}.dcrkt,1);
        for iapp = 1:n_app
            % dcrkt
            tmp = app_f{icond,ilaser}.dcrkt(iapp,:);
            tmp_nan = isnan(tmp);
            tmp(tmp_nan) = [];
            tmp_len = length(tmp);

            if tmp_len > 5
                tmp_x = 1:tmp_len;
                tmp_y = tmp;
                tmp_interp = 1:((tmp_len-1)/smp):tmp_len;
    
                vq = interpn(tmp_x,tmp_y,tmp_interp,'spline');
    
                interp_f{icond,ilaser}.dcrkt = [interp_f{icond,ilaser}.dcrkt; vq];
            end
        end

        n_app = size(app_f{icond,ilaser}.cspd,1);
        for iapp = 1:n_app
            % cspd
            tmp = app_f{icond,ilaser}.cspd(iapp,:);
            tmp_nan = isnan(tmp);
            tmp(tmp_nan) = [];
            tmp_len = length(tmp);

            if tmp_len > 5
                tmp_x = 1:tmp_len;
                tmp_y = tmp;
                tmp_interp = 1:((tmp_len-1)/smp):tmp_len;
    
                vq = interpn(tmp_x,tmp_y,tmp_interp,'spline');
    
                interp_f{icond,ilaser}.cspd = [interp_f{icond,ilaser}.cspd; vq];
            end
        end
    end
end

% plot:

cindex = [0 0.75 0.75; 0 0 0];
for iparam = 1:4

    f = figure(3777+iparam); clf
    f.Position = [1957 199 1031 1469];

    ylimcmp = [];
    for icond = 1:8
        for ilaser = 1:2

            if iparam == 1
                tmp_s = interp_s{icond,ilaser}.spd;
                tmp_f = interp_f{icond,ilaser}.spd;
            elseif iparam == 2
                tmp_s = interp_s{icond,ilaser}.az;
                tmp_f = interp_f{icond,ilaser}.az;
            elseif iparam == 3
                tmp_s = interp_s{icond,ilaser}.dcrkt;
                tmp_f = interp_f{icond,ilaser}.dcrkt;
            elseif iparam == 4
                tmp_s = interp_s{icond,ilaser}.cspd;
                tmp_f = interp_f{icond,ilaser}.cspd;
            end

            s_x = 1:101;
            s_y = mean(tmp_s);
            s_e = std(tmp_s)/sqrt(size(tmp_s,1));
            subplot(8,2,(icond-1)*2+ilaser); cla
            S1 = shadedErrorBar(s_x,s_y,s_e); hold on
            S1.patch.FaceColor = cindex(1,:);

            s_x = 1:101;
            s_y = mean(tmp_f);
            s_e = std(tmp_f)/sqrt(size(tmp_f,1));
            subplot(8,2,(icond-1)*2+ilaser);
            S2 = shadedErrorBar(s_x,s_y,s_e); hold on
            S2.patch.FaceColor = cindex(2,:);

            if ilaser == 1
                title([cond_list{icond} ', Loff'])
            elseif ilaser == 2
                title([cond_list{icond} ', Lon'])
            end

            % if icond == 1
            %     legend('success','fail')
            % end
            ax1 = gca;
            ylimcmp = [ylimcmp; ax1.YLim(1) ax1.YLim(2)];
        end
    end

    for icond = 1:8
        for ilaser = 1:2
            ylimmin = min(ylimcmp(:,1));
            ylimmax = max(ylimcmp(:,2));

            subplot(8,2,(icond-1)*2+ilaser);

            xlim([0 101])
            ylim([ylimmin ylimmax])
        end
    end

    if iparam == 1
        sav_title = ['norm_spd.png'];
    elseif iparam == 2
        sav_title = ['norm_az.png'];
    elseif iparam == 3
        sav_title = ['norm_dcrkt.png'];
    elseif iparam == 4
        sav_title = ['norm_cspd.png'];
    end


    F = getframe(f);
    imwrite(F.cdata, sav_title, 'png')
    close
end