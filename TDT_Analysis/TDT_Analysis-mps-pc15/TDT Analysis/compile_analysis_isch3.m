%function compile_analysis_movement
clear all

N=0;

% CSD experiments: 90625, 90728, 90806-3, 

% CONTROL EXPERIMENTS
afiles(1,:) = {'C:\data\isch\90622 analysis\analysis.mat' ...
    'C:\data\isch\90623 analysis\analysis_nomerge.mat' ...
    'C:\data\isch\90625 analysis\analysis_nomerge.mat' ... % cluster...123_2
    'C:\data\isch\90728 analysis2\analysis45.mat' ... 
    'C:\data\isch\90728 analysis2\analysis1112.mat' ...
    };

% HI EXPERIMENTS
afiles(2,:) = {'C:\data\isch\90518 analysis\analysis2_isch111210.mat' ... % injury = mild (?)
    'C:\data\isch\90604_0907-1 analysis\analysis.mat' ... % injury = unknown (using 1)
    'C:\data\isch\90605_0907-11 analysis\analysis.mat' ... % injury = moderate     % cluster...isch23
    'C:\data\isch\90617_0907-5 analysis\analysis.mat' ... % injury = moderate (no hist?)
    'C:\data\isch\90806-3 analysis\analysis.mat' ... % injury = moderate
    };

for (fileidx=1:2) % 1 = CONTROL, 2 = INJURED
    for i = 1:length(afiles(fileidx,:))
        load(afiles{fileidx,i});
        afiles(fileidx,i);

        n_units = length(L_ratio);
        cellrange=N+1:N+n_units;
        N=N+n_units;

        ch_num(cellrange)=cells(:,1);
        cl_num(cellrange)=cells(:,2);
        for j = cellrange
            fname{j} = afiles(fileidx,i);
        end

        if (fileidx==1) % CONTROLS
            injured_allcells(cellrange)=0;
        else
            injured_allcells(cellrange)=injured;
        end

        drift_peak_allcells(cellrange)=drift_peak;
        drift_spont_allcells(cellrange)=drift_spont;
        drift_OSI_allcells(cellrange)=drift_OSI;
        drift_thetawidth_allcells(cellrange)=drift_thetawidth;
            
        r_f0_allcells(cellrange)=r_f0;
        r_f1_allcells(cellrange)=r_f1;
        r_f1_phase_allcells(cellrange)=r_f1_phase;
        r_f1dividef0_allcells(cellrange)=r_f1dividef0;

        trough2peak_allcells(cellrange)=trough2peak;

        wv_allcells(cellrange,:)=wv';
    end
end

% ALIGN WAVEFORMS

wv=wv_allcells;
for x=1:size(wv,1)
    %         if inh(i)==1;
    %         colorstr='b';
    %     elseif inh(i)==0;
    %         colorstr='g';
    %     else colorstr='w';
    %         end
    %     if ~used(i)
    %         colorstr='w';
    %     end
    colorstr = 'k';
    [y j] = min(wv(x,:));
    j;
    x;
    if j ==5
        wvshift(x,:) = wv(x,1:17);
%        plot(wv(x,1:17),colorstr);
    elseif j==6
        wvshift(x,:) = wv(x,2:18);
%        plot(wv(x,2:18),colorstr);
    elseif j==7
        wvshift(x,:) = wv(x,3:19);
%        plot(wv(x,3:19),colorstr);
    else
        wvshift(x,:) = wv(x,1:17);
        used(x)=0;
    end
end

used = ones(N,1)';
used(wvshift(:,14)<0.3540)=0; % filter out a few units with non-spike waveforms
used(23)=0;
used(73)=0;

% PLOT INH/EXC CELLS
new_endheight = wvshift(:,17)-wvshift(:,14);
[new_max t2peak]= max(wvshift,[],2);

inh = trough2peak_allcells<9.5;
used(inh)=0;

% hold on
% for i=1:size(wv,1)
% %     if inh(i)==1;
% %         colorstr='b';
% %     elseif inh(i)==0;
% %         colorstr='g';
% %     end
%     if injured_allcells(i)>0
%             colorstr='r';
%     end
%     
%     if used(1,i)==1;
%         plot(wvshift(i,1:17),colorstr);
%     end
% end

% OUTPUT VARIABLES
output=1;
if (output==1)
    usedcells = find(used == 1);   
    
    % map [-pi, pi] to [0, 2*pi] -- r_f1_rephase
    for i = 1:length(r_f1_phase_allcells)
        if (r_f1_phase_allcells(i) > 0) 
            r_f1_rephase_allcells(i)=r_f1_phase_allcells(i)-2*pi
        else
            r_f1_rephase_allcells(i)=r_f1_phase_allcells(i)
        end
    end
   
    meanphase=angle(mean(cos(r_f1_rephase_allcells(usedcells)+sqrt(-1)*sin(r_f1_rephase_allcells(usedcells)))));
    phaserange=30*pi/180; % +/- 30 degrees

    % CONTROLS
    con=find(injured_allcells==0 & used==1);
    con_resp=find(drift_peak_allcells>2 & injured_allcells==0 & used==1);
    con_phase=find(r_f1_rephase_allcells>meanphase-phaserange & r_f1_rephase_allcells<meanphase+phaserange & injured_allcells==0 & used==1);
    con_nophase=find((r_f1_rephase_allcells<meanphase-phaserange | r_f1_rephase_allcells>meanphase+phaserange) & injured_allcells==0 & used==1);

    con_orient = find(drift_peak_allcells>2 & drift_OSI_allcells>0.5 & injured_allcells==0 & used==1);
  
    numof_con=length(con);
    numof_resp_con=length(con_resp);
    numof_phase_con=length(con_phase);
    numof_orient_con=length(con_orient);

    percent_resp_con=(numof_resp_con/numof_con)*100;
    percent_phase_con=(numof_phase_con/numof_con)*100;
    percent_orient_con=(numof_orient_con/numof_resp_con)*100;

    percent_ch1_con=length(find(ch_num(con)==1))/numof_con*100;
    percent_ch5_con=length(find(ch_num(con)==5))/numof_con*100;
    percent_ch9_con=length(find(ch_num(con)==9))/numof_con*100;
    percent_ch13_con=length(find(ch_num(con)==13))/numof_con*100;
    
    mean_drift_spont_con=mean(drift_spont_allcells(con));
    mean_drift_spont_con_SEM=std(drift_spont_allcells(con))/sqrt(length(drift_spont_allcells(con)));
    mean_drift_OSI_con=mean(drift_OSI_allcells(con_resp));
    mean_drift_OSI_con_SEM=std(drift_OSI_allcells(con_resp))/sqrt(length(drift_OSI_allcells(con_resp)));
    mean_drift_thetawidth_con=mean(drift_thetawidth_allcells(con_orient));
    mean_drift_thetawidth_con_SEM=std(drift_thetawidth_allcells(con_orient))/sqrt(length(drift_thetawidth_allcells(con_orient)));
    mean_f1all_con=mean(r_f1_allcells(con));
    mean_f1all_con_SEM=std(r_f1_allcells(con))/sqrt(length(r_f1_allcells(con)));
    mean_f1phase_con=mean(r_f1_allcells(con_phase));
    mean_f1phase_con_SEM=std(r_f1_allcells(con_phase))/sqrt(length(r_f1_allcells(con_phase)));
    
    median_drift_spont_con=median(drift_spont_allcells(con));
    median_drift_spont_con_SEM=std(drift_spont_allcells(con))/sqrt(length(drift_spont_allcells(con)));
    median_drift_OSI_con=median(drift_OSI_allcells(con_resp));
    median_drift_OSI_con_SEM=std(drift_OSI_allcells(con_resp))/sqrt(length(drift_OSI_allcells(con_resp)));
    median_drift_thetawidth_con=median(drift_thetawidth_allcells(con_orient));
    median_drift_thetawidth_con_SEM=std(drift_thetawidth_allcells(con_orient))/sqrt(length(drift_thetawidth_allcells(con_orient)));
    median_f1all_con=median(r_f1_allcells(con));
    median_f1all_con_SEM=std(r_f1_allcells(con))/sqrt(length(r_f1_allcells(con)));
    median_f1phase_con=median(r_f1_allcells(con_phase));
    median_f1phase_con_SEM=std(r_f1_allcells(con_phase))/sqrt(length(r_f1_allcells(con_phase)));

    % INJURED
    hi=find(injured_allcells>0 & used==1);
    hi_resp=find(drift_peak_allcells>2 & injured_allcells>0 & used==1);
    hi_phase=find(r_f1_rephase_allcells>meanphase-phaserange & r_f1_rephase_allcells<meanphase+phaserange & injured_allcells>0 & used==1);   
    hi_nophase=find((r_f1_rephase_allcells<meanphase-phaserange | r_f1_rephase_allcells>meanphase+phaserange) & injured_allcells>0 & used==1);   
    hi_orient = find(drift_peak_allcells>2 & drift_OSI_allcells>0.5 & injured_allcells>0 & used==1);
  
    numof_hi=length(hi);
    numof_resp_hi=length(hi_resp);
    numof_phase_hi=length(hi_phase);
    numof_orient_hi=length(hi_orient);

    percent_resp_hi=(numof_resp_hi/numof_hi)*100;
    percent_phase_hi=(numof_phase_hi/numof_hi)*100;
    percent_orient_hi=(numof_orient_hi/numof_resp_hi)*100;

    percent_ch1_hi=length(find(ch_num(hi)==1))/numof_hi*100;
    percent_ch5_hi=length(find(ch_num(hi)==5))/numof_hi*100;
    percent_ch9_hi=length(find(ch_num(hi)==9))/numof_hi*100;
    percent_ch13_hi=length(find(ch_num(hi)==13))/numof_hi*100;
    
    mean_drift_spont_hi=mean(drift_spont_allcells(hi));
    mean_drift_spont_hi_SEM=std(drift_spont_allcells(hi))/sqrt(length(drift_spont_allcells(hi)));
    mean_drift_OSI_hi=mean(drift_OSI_allcells(hi_resp));
    mean_drift_OSI_hi_SEM=std(drift_OSI_allcells(hi_resp))/sqrt(length(drift_OSI_allcells(hi_resp)));
    mean_drift_thetawidth_hi=mean(drift_thetawidth_allcells(hi_orient));
    mean_drift_thetawidth_hi_SEM=std(drift_thetawidth_allcells(hi_orient))/sqrt(length(drift_thetawidth_allcells(hi_orient)));
    mean_f1all_hi=mean(r_f1_allcells(hi));
    mean_f1all_hi_SEM=std(r_f1_allcells(hi))/sqrt(length(r_f1_allcells(hi)));
    mean_f1phase_hi=mean(r_f1_allcells(hi_phase));
    mean_f1phase_hi_SEM=std(r_f1_allcells(hi_phase))/sqrt(length(r_f1_allcells(hi_phase)));

    median_drift_spont_hi=median(drift_spont_allcells(hi));
    median_drift_spont_hi_SEM=std(drift_spont_allcells(hi))/sqrt(length(drift_spont_allcells(hi)));
    median_drift_OSI_hi=median(drift_OSI_allcells(hi_resp));
    median_drift_OSI_hi_SEM=std(drift_OSI_allcells(hi_resp))/sqrt(length(drift_OSI_allcells(hi_resp)));
    median_drift_thetawidth_hi=median(drift_thetawidth_allcells(hi_orient));
    median_drift_thetawidth_hi_SEM=std(drift_thetawidth_allcells(hi_orient))/sqrt(length(drift_thetawidth_allcells(hi_orient)));
    median_f1all_hi=median(r_f1_allcells(hi));
    median_f1all_hi_SEM=std(r_f1_allcells(hi))/sqrt(length(r_f1_allcells(hi)));
    median_f1phase_hi=median(r_f1_allcells(hi_phase));
    median_f1phase_hi_SEM=std(r_f1_allcells(hi_phase))/sqrt(length(r_f1_allcells(hi_phase)));

    countvars = [ {'numof'} numof_con numof_hi; ...
        {'numof_resp'} numof_resp_con numof_resp_hi; ...
        {'numof_phase'} numof_phase_con numof_phase_hi; ...
        {'numof_orient'} numof_orient_con numof_orient_hi; ...
        {'percent_resp'} percent_resp_con percent_resp_hi; ...
        {'percent_phase'} percent_phase_con percent_phase_hi; ...
        {'percent_orient'} percent_orient_con percent_orient_hi; ...
        {'percent_ch1'} percent_ch1_con percent_ch1_hi; ...
        {'percent_ch5'} percent_ch5_con percent_ch5_hi; ...
        {'percent_ch9'} percent_ch9_con percent_ch9_hi; ...
        {'percent_ch13'} percent_ch13_con percent_ch13_hi ]
    meanvars = [ {'mean_drift_spont'} mean_drift_spont_con mean_drift_spont_hi mean_drift_spont_con_SEM mean_drift_spont_hi_SEM; ...
        {'mean_drift_OSI'} mean_drift_OSI_con mean_drift_OSI_hi mean_drift_OSI_con_SEM mean_drift_OSI_hi_SEM; ...
        {'mean_drift_thetawidth'} mean_drift_thetawidth_con mean_drift_thetawidth_hi mean_drift_thetawidth_con_SEM mean_drift_thetawidth_hi_SEM; ...
        {'mean_f1all'} mean_f1all_con mean_f1all_hi mean_f1all_con_SEM mean_f1all_hi_SEM; ...
        {'mean_f1phase'} mean_f1phase_con mean_f1phase_hi mean_f1phase_con_SEM mean_f1phase_hi_SEM ]
    medianvars = [ {'median_drift_spont'} median_drift_spont_con median_drift_spont_hi median_drift_spont_con_SEM median_drift_spont_hi_SEM; ...
        {'median_drift_OSI'} median_drift_OSI_con median_drift_OSI_hi median_drift_OSI_con_SEM median_drift_OSI_hi_SEM; ...
        {'median_drift_thetawidth'} median_drift_thetawidth_con median_drift_thetawidth_hi median_drift_thetawidth_con_SEM median_drift_thetawidth_hi_SEM; ...
        {'median_f1all'} median_f1all_con median_f1all_hi median_f1all_con_SEM median_f1all_hi_SEM; ...
        {'median_f1phase'} median_f1phase_con median_f1phase_hi median_f1phase_con_SEM median_f1phase_hi_SEM ]

    % STATISTICAL COMPARISON
    'spontaneous firing rates'
    [P,H] = ranksum(drift_spont_allcells(con), drift_spont_allcells(hi))
    [H,P] = ttest2(drift_spont_allcells(con), drift_spont_allcells(hi))
    
    'orientation selectivity index (OSI)'
    [P,H] = ranksum(drift_OSI_allcells(con_resp), drift_OSI_allcells(hi_resp))
    [H,P] = ttest2(drift_OSI_allcells(con_resp), drift_OSI_allcells(hi_resp))
    
    'drift thetawidth'
    [P,H] = ranksum(drift_thetawidth_allcells(con_orient), drift_thetawidth_allcells(hi_orient))
    [H,P] = ttest2(drift_thetawidth_allcells(con_orient), drift_thetawidth_allcells(hi_orient))
    
    'f1 noise movie response (all cells)'
    [P,H] = ranksum(r_f1_allcells(con), r_f1_allcells(hi))
    [H,P] = ttest2(r_f1_allcells(con), r_f1_allcells(hi))
    
    'f1 noise movie response (phase matched cells)'
    [P,H]=ranksum(r_f1_allcells(con_phase),r_f1_allcells(hi_phase))
    [H,P]=ttest2(r_f1_allcells(con_phase),r_f1_allcells(hi_phase))

    % average waveforms
    wvshift_con=mean(wvshift(con,1:17));
    wvshift_hi=mean(wvshift(hi,1:17));

    meanphase_con=angle(mean(cos(r_f1_rephase_allcells(con_phase)+sqrt(-1)*sin(r_f1_rephase_allcells(con_phase)))));
    meanphase_hi=angle(mean(cos(r_f1_rephase_allcells(hi_phase)+sqrt(-1)*sin(r_f1_rephase_allcells(hi_phase)))));

    meanamp_nophase_con=mean(r_f1_allcells(con_nophase))
    meanamp_nophase_hi=mean(r_f1_allcells(hi_nophase))

    markerphase=[0 meanphase-phaserange 0 meanphase+phaserange];
    markeramp=[0 5 0 5];

    polar(markerphase, markeramp, 'g')
    hold on
    polar([0 meanphase_con], [0 mean(r_f1_allcells(con_phase))], 'k') % con line
    polar([0 meanphase_hi], [0 mean(r_f1_allcells(hi_phase))], 'r') % hi line
    polar(meanphase_con, mean(r_f1_allcells(con_phase)), '*k') % con point
    polar(meanphase_hi, mean(r_f1_allcells(hi_phase)), '*r') % hi point
    
    % all points
    polar(r_f1_rephase_allcells(con), r_f1_allcells(con), '.g');
    polar(r_f1_rephase_allcells(hi), r_f1_allcells(hi), '.g');
    polar(r_f1_rephase_allcells(con_phase), r_f1_allcells(con_phase), '.k');
    polar(r_f1_rephase_allcells(hi_phase), r_f1_allcells(hi_phase), '.r');
end

% XLS EXPORT
xlsoutput=0;
if (xlsoutput==1)    
    xlsfile='c:\data\isch\analysis_allcells.xls'

    % extract cell contents of fname into fname_out
    for (i=1:size(fname,2))
        fname_out(i)=fname{i};
    end
    xlswrite(xlsfile, ch_num(con_resp)', 'con_resp cells', 'A1')
    xlswrite(xlsfile, cl_num(con_resp)', 'con_resp cells', 'B1')
    xlswrite(xlsfile, fname_out(con_resp)', 'con_resp cells', 'C1')
    xlswrite(xlsfile, ch_num(hi_resp)', 'hi_resp cells', 'A1')
    xlswrite(xlsfile, cl_num(hi_resp)', 'hi_resp cells', 'B1')
    xlswrite(xlsfile, fname_out(hi_resp)', 'hi_resp cells', 'C1')
    xlswrite(xlsfile, drift_spont_allcells(con)', 'drift_spont_allcells', 'A1')
    xlswrite(xlsfile, drift_spont_allcells(hi)', 'drift_spont_allcells', 'B1')
    xlswrite(xlsfile, drift_OSI_allcells(con_resp)', 'drift_OSI_resp', 'A1')
    xlswrite(xlsfile, drift_OSI_allcells(hi_resp)', 'drift_OSI_resp', 'B1')
    xlswrite(xlsfile, drift_thetawidth_allcells(con_orient)', 'drift_thetawidth_orient', 'A1')
    xlswrite(xlsfile, drift_thetawidth_allcells(hi_orient)', 'drift_thetawidth_orient', 'B1')
    xlswrite(xlsfile, r_f1_allcells(con)', 'f1_allcells', 'A1')
    xlswrite(xlsfile, r_f1_allcells(hi)', 'f1_allcells', 'B1')
    xlswrite(xlsfile, r_f1_allcells(con_resp)', 'f1_resp', 'A1')
    xlswrite(xlsfile, r_f1_allcells(hi_resp)', 'f1_resp', 'B1')
    xlswrite(xlsfile, r_f1_allcells(con_phase)', 'f1_phase', 'A1')
    xlswrite(xlsfile, r_f1_allcells(hi_phase)', 'f1_phase', 'B1')
    xlswrite(xlsfile, countvars, 'vars', 'A2')
    xlswrite(xlsfile, meanvars, 'vars', 'D2')   
    xlswrite(xlsfile, medianvars, 'vars', 'D10')


end



savefile='c:\data\isch\analysis_allcells'
save(savefile, 'injured_allcells', ...
    'r_f0_allcells','r_f1_allcells','r_f1_phase_allcells','r_f1dividef0_allcells', ...
    'drift_peak_allcells','drift_spont_allcells','drift_OSI_allcells' ...
    );

%end % main function