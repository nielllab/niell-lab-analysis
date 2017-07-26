ON=1; OFF=0;
%pathname = 'F:Prey_tracks\';%miyazaki
pathname='X:\';% miyazaki path
%pathname='E:\Prey_capture_data\'; %eyre path
if ismac
    pathname = '/Users/crisniell/Dropbox/Prey capture/plexi data/'
end


n=0;

%%%% prey 1
%%%group: 1=Plexiglass 1 side; 2=plexiglass 2 side 


n=n+1;
files(n).subj = 'Prey12_PT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_5_31_16\Prey_12_PT_trial_1_LS_NIR\Tracks\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=6;%blkBgd but included in ear plug group white background
files(n).sex=2;
files(n).FrameS=3;%frame where cricket is first available
files(n).FrameEnd=3005%sample,only 4188 points included in track data, 7200 long total video
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_whitebgd_4_5_16\Prey_12_PT_trial_1_LS_NIR.avi';
files(n).body=0;
files(n).EP=1;



n=n+1;
files(n).subj = 'Prey7_TT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_whitebgd_4_5_16\Prey7_TTearplug_trial_4_leftside_NIR\Tracks\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3000;
files(n).FrameFlash=''
files(n).Moviefile='Prey7_TTearplug_trial_4_leftside_NIR.avi'
files(n).EP=1;

n=n+1;
files(n).subj = 'Prey8-RT';
files(n).lighting = ON;
files(n).Tnum=2;
files(n).trackpts = 'Plexiglass_3_18_16\Prey8_RT_PG_FS_choice_trial4\Tracks\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=2;
files(n).FrameS=1;
files(n).FrameEnd=3000;
files(n).FrameFlash=''
files(n).Moviefile='Plexiglass_3_18_16\Prey8_RT_PG_FS_choice_trial4.avi'
files(n).EP=1;

n=n+1;
files(n).subj = 'Prey7-NT';
files(n).lighting = OFF;
files(n).Tnum=2;
files(n).trackpts = 'Plexiglass_3_18_16\Prey7_NT_PG_NIR_choice_trial5\tracks\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3000;
files(n).FrameFlash=''
files(n).Moviefile='Plexiglass_3_18_16\Prey7_NT_PG_NIR_choice_trial5.avi'
files(n).EP=1;

n=n+1;
files(n).subj = 'Prey7-NT';
files(n).lighting = ON;
files(n).Tnum=2;
files(n).trackpts = 'Plexiglass_3_18_16\Prey7_NT_FS_choice_trial4\Tracks\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3000;
files(n).FrameFlash=''
files(n).Moviefile='Plexiglass_3_18_16\Prey7_NT_FS_choice_trial4.avi'
files(n).EP=1;

n=n+1;
files(n).subj = 'Prey10-LT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_3_18_16\Prey10_LT_PG_FS_LS_trial4\Tracks\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=2;
files(n).FrameS=1;
files(n).FrameEnd=3000;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_3_18_16\Prey10_LT_PG_FS_LS_trial4.avi';
files(n).EP=1;

n=n+1;
files(n).subj = 'Prey11-RT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_3_18_16\Prey11_RT_PG_FS_trial1_RS\Tracks\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3000;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_3_18_16\Prey11_RT_PG_FS_trial1_RS.avi';
files(n).EP=1;

n=n+1;
files(n).subj = 'Prey11-RT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_3_18_16\Prey11_RT_PG_FS_trial3_LS\Tracks\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3000;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_3_18_16\Prey11_RT_PG_FS_trial3_LS.avi';
files(n).EP=1;

n=n+1;
files(n).subj = 'Prey11-RT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_3_18_16\Prey11_RT_PG_NIR_trial2_LS\Tracks\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3000;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_3_18_16\Prey11_RT_PG_NIR_trial2_LS.avi';
files(n).EP=1;

n=n+1;
files(n).subj = 'Prey12-NT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey_12_NT_trial_1_FS_LS\Tracks\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=2;
files(n).FrameS=1;
files(n).FrameEnd=3000;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_whitebgd_4_15_16\Prey_12_NT_trial_1_FS_LS.avi';
files(n).EP=1;

n=n+1;
files(n).subj = 'Prey12_PT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey_12_PT_trial_4_FS_LS\Tracks\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=2;
files(n).FrameS=1;
files(n).FrameEnd=3000;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_whitebgd_4_15_16\Prey_12_PT_trial_4_FS_LS.avi';
files(n).EP=1;

n=n+1;
files(n).subj = 'Prey13_LT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey_13_LT_trial_1_NIR_RS\Tracks\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=2;
files(n).FrameS=3;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_whitebgd_4_15_16\BW_Prey_13_LT_trial_1_NIR_RS.avi';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey13_LT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey_13_LT_trial_2_FS_RS\Tracks\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=2;
files(n).FrameS=3;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_whitebgd_4_15_16\Prey_13_LT_trial_2_FS_RS.avi';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey13_LT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey_13_LT_trial_3_NIR_LS\Tracks\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=2;
files(n).FrameS=3;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_whitebgd_4_15_16\Prey_13_LT_trial_3_NIR_LS.avi';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey13_LT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey_13_LT_trial_4_FS_LS\Tracks\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=2;
files(n).FrameS=3;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_whitebgd_4_15_16\Prey_13_LT_trial_4_FS_LS.avi';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey14_NT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey_14_NT_trial_1_NIR_LS\Tracks\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=3;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_whitebgd_4_15_16\Prey_14_NT_trial_1_NIR_LS.avi';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey14_NT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey_14_NT_trial_2_FS_LS\Tracks\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=3;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_whitebgd_4_15_16\Prey_14_NT_trial_2_FS_LS.avi';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey14_NT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey_14_NT_trial_4_FS_RS\Tracks\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=3;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_whitebgd_4_15_16\Prey_14_NT_trial_4_FS_RS.avi';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey12_PT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_5_31_16\Prey12_PT_trial3_NIR\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_5_31_16\Prey_12_PT_trial_3_RS_NIR.avi';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey13_NT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_5_31_16\Prey13_NT_trial4_NIR\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_5_31_16\Prey_13_NT_trial_4_RS_NIR.avi';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey13_PT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_5_31_16\Prey13_PT_trial4_NIR\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_5_31_16\Prey_13_PT_trial_4_LS_NIR.avi';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey14_PT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_5_31_16\Prey_14_PT_trial1_NIR\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_5_31_16\Prey_14_PT_trial_1_RS_NIR.avi';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey14_PT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_5_31_16\Prey14_PT_trial3_NIR\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_5_31_16\Prey_14_PT_trial_3_LS_NIR.avi';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey12_PT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_5_31_16\Prey12_PT_trail1_NIR_2\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_5_31_16\Prey_14_PT_trial_3_LS_NIR.avi';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey7_TT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_whitebgd_4_5_16\Prey7_TT_trial4_NIR\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey4_LT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Old_plexi_trials\Prey4_LT_trial2_NIR\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_whitebgd_4_5_16\Prey7_TTearplug_trial_4_leftside_NIR.avi';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey4_LT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Old_plexi_trials\Prey4_LT_trial4_NIR\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey5_LT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_5_31_16\Prey5_LT_trial2_NIR\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_whitebgd_4_5_16\Prey7_TTearplug_trial_4_leftside_NIR.avi';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey5_RT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_5_31_16\Prey5_LT_NIR_trial6\DLTdv5_data_xypts.csv';%trial6 Prey5RT error in folders?
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='Plexiglass_whitebgd_4_5_16\Prey7_TTearplug_trial_4_leftside_NIR.avi';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey5_RT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_5_31_16\Prey5_RT_trial3\DLTdv5_data_xypts.csv';%trial6 Prey5RT error in folders?
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd='';
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey4_LT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Old_plexi_trials\Prey4_LT_trial3_NIR\DLTdv5_data_xypts.csv';%trial6 Prey5RT error in folders?
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd='';
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey7_TT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_whitebgd_4_5_16\Prey7_TT_NIR_trial1\DLTdv5_data_xypts.csv';%trial6 Prey5RT error in folders?
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd='';
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey7_TT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_whitebgd_4_5_16\Prey7_TT_NIR_trial2\DLTdv5_data_xypts.csv';%trial6 Prey5RT error in folders?
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd='';
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey7_NT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_whitebgd_4_5_16\Prey7_NT_NIR_trial2\DLTdv5_data_xypts.csv';%trial6 Prey5RT error in folders?
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd='';
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey13_NT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey12_NT_trial3_FS\DLTdv5_data_xypts.csv';%trial6 Prey5RT error in folders?
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd='';
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey12_PT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_5_31_16\Light\Prey12_PT_trial4_FS\DLTdv5_data_xypts.csv';%trial6 Prey5RT error in folders?
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd='';
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey15_LT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_5_31_16\Light\Prey5_LT_trial1\DLTdv5_data_xypts.csv';%trial6 Prey5RT error in folders?
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd='';
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey15_LT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_5_31_16\Light\Prey5_LT_trial6\DLTdv5_data_xypts.csv';%trial6 Prey5RT error in folders?
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd='';
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey14_NT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey14_NT_trial3_FS\DLTdv5_data_xypts.csv';%trial6 Prey5RT error in folders?
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd='';
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey13_LT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey13_LT_trial6_FS\DLTdv5_data_xypts.csv';%trial6 Prey5RT error in folders?
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd='';
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey13_LT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey13_LT_trial7\DLTdv5_data_xypts.csv';%trial6 Prey5RT error in folders?
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd='';
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey14_NT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey14_NT_trial2_FS_3\DLTdv5_data_xypts.csv';%trial6 Prey5RT error in folders?
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd='';
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey14_NT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey14_NT_trial4_T2\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey14_PT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_5_31_16\Prey14_PT_trial3_T2\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey12_PT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexiglass_5_31_16\Prey12_PT_trial1_NIR_T3\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey12_NT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexi_6_20_16\Prey12_NT_FS_LS_trial3\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey12_NT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexi_6_20_16\Prey12_NT_FS_RS_trial1\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey12_NT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexi_6_20_16\Prey12_NT_NIR_LS_Trial4\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey12_NT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexi_6_20_16\Prey12_NT_NIR_RS_trial2\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey12_pT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexi_6_20_16\Prey12_pT_FS_RS_trial2\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey12_pT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexi_6_20_16\Prey12_pT_NIR_LS_trial3\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey12_pT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexi_6_20_16\Prey12_pT_NIR_RS_trial1\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey16_pT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexi_6_20_16\Prey16_NT_FS_LS_Trial1\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey16_NT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexi_6_20_16\Prey16_NT_FS_LS_trial3\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey16_NT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexi_6_20_16\Prey16_NT_NIR_LS_trial2\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey16_pT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexi_6_20_16\Prey16_pT_FS_LS_trial1\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey16_pT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexi_6_20_16\Prey16_pT_FS_RS_trial3\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey16_pT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexi_6_20_16\Prey16_pT_NIR_LS_trial4\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey17_NT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexi_6_20_16\Prey17_NT_FS_LS_trial2b\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey17_NT';
files(n).lighting = ON;
files(n).Tnum=1;
files(n).trackpts = 'Plexi_6_20_16\Prey17_NT_FS_RS_Trial4\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey17_NT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexi_6_20_16\Prey17_NT_NIR_LS_trial1B\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey17_NT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexi_6_20_16\Prey17_NT_NIR_RS_trial3\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;


n=n+1;
files(n).subj = 'Prey17_NT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexi_6_20_16\Prey17_pT_NIR_LS_trial1\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;

n=n+1;
files(n).subj = 'Prey17_NT';
files(n).lighting = OFF;
files(n).Tnum=1;
files(n).trackpts = 'Plexi_6_20_16\Prey17_pT_NIR_RS_trial3\DLTdv5_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=1;
files(n).FrameEnd=3600;
files(n).FrameFlash='';
files(n).Moviefile='';
files(n).EP=0;


