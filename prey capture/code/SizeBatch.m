ON=1; OFF=0;
pathname = 'X:\Preycapture\';%miyazaki path to prey capture folder stored on Eyre 
n=0;


%%%group_contrast: 1=100; 2=50 ; 3=25; 4=12.5; 5=6.25;

%%%group_size:1=1x, 2=2x,3=4x, 4=.5x, 5=.25x, 6=.125x, 7=n1x, 8=n2x, 9=n4x, 10=.5x,
%%%11=.25x, 12=.125x

%%%female=1 and male=2;
%%%contrast:100,50, 25, 12.5, and 6.25, percent transmittance
%%%size:100, 50, 25, 12.5, 6.25, scaled down size along all dimensions of
%%%ellipse (naturalistic stimuli)

%%
% n=n+1;
% files(n).subj = 'Prey12_PT';
% files(n).lighting = OFF;
% files(n).trackpts = 'Plexiglass_5_31_16\Prey_12_PT_trial_1_LS_NIR\Tracks\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;%blkBgd but included in ear plug group white background
% files(n).Moviefile='Plexiglass_whitebgd_4_5_16\Prey_12_PT_trial_1_LS_NIR.avi';
% 

% n=n+1;
% files(n).subj = 'Prey7_TT';
% files(n).lighting = OFF;
% files(n).trackpts = 'Plexiglass_whitebgd_4_5_16\Prey7_TTearplug_trial_4_leftside_NIR\Tracks\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Prey7_TTearplug_trial_4_leftside_NIR.avi'
% 
% 
% n=n+1;
% files(n).subj = 'Prey8-RT';
% files(n).lighting = ON;
% files(n).trackpts = 'Plexiglass_3_18_16\Prey8_RT_PG_FS_choice_trial4\Tracks\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_3_18_16\Prey8_RT_PG_FS_choice_trial4.avi'
% 
% 
% n=n+1;
% files(n).subj = 'Prey7-NT';
% files(n).lighting = OFF;
% files(n).trackpts = 'Plexiglass_3_18_16\Prey7_NT_PG_NIR_choice_trial5\tracks\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_3_18_16\Prey7_NT_PG_NIR_choice_trial5.avi'
% 
% 
% n=n+1;
% files(n).subj = 'Prey7-NT';
% files(n).lighting = ON;
% files(n).trackpts = 'Plexiglass_3_18_16\Prey7_NT_FS_choice_trial4\Tracks\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_3_18_16\Prey7_NT_FS_choice_trial4.avi'
% 
% 
% n=n+1;
% files(n).subj = 'Prey10-LT';
% files(n).lighting = ON;
% files(n).Tnum=1;
% files(n).trackpts = 'Plexiglass_3_18_16\Prey10_LT_PG_FS_LS_trial4\Tracks\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_3_18_16\Prey10_LT_PG_FS_LS_trial4.avi';
% 
% 
% n=n+1;
% files(n).subj = 'Prey11-RT';
% files(n).lighting = ON;
% files(n).trackpts = 'Plexiglass_3_18_16\Prey11_RT_PG_FS_trial1_RS\Tracks\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_3_18_16\Prey11_RT_PG_FS_trial1_RS.avi';
% 
% 
% n=n+1;
% files(n).subj = 'Prey11-RT';
% files(n).lighting = ON;
% files(n).trackpts = 'Plexiglass_3_18_16\Prey11_RT_PG_FS_trial3_LS\Tracks\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_3_18_16\Prey11_RT_PG_FS_trial3_LS.avi';
% 
% 
% n=n+1;
% files(n).subj = 'Prey11-RT';
% files(n).lighting = OFF;
% files(n).trackpts = 'Plexiglass_3_18_16\Prey11_RT_PG_NIR_trial2_LS\Tracks\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_3_18_16\Prey11_RT_PG_NIR_trial2_LS.avi';
% 
% 
% n=n+1;
% files(n).subj = 'Prey12-NT';
% files(n).lighting = ON;
% files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey_12_NT_trial_1_FS_LS\Tracks\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_whitebgd_4_15_16\Prey_12_NT_trial_1_FS_LS.avi';
% 
% 
% n=n+1;
% files(n).subj = 'Prey12_PT';
% files(n).lighting = ON;
% files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey_12_PT_trial_4_FS_LS\Tracks\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_whitebgd_4_15_16\Prey_12_PT_trial_4_FS_LS.avi';
% 
% 
% n=n+1;
% files(n).subj = 'Prey13_LT';
% files(n).lighting = OFF;
% files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey_13_LT_trial_1_NIR_RS\Tracks\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_whitebgd_4_15_16\BW_Prey_13_LT_trial_1_NIR_RS.avi';
% 
% 
% n=n+1;
% files(n).subj = 'Prey13_LT';
% files(n).lighting = ON;
% files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey_13_LT_trial_2_FS_RS\Tracks\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_whitebgd_4_15_16\Prey_13_LT_trial_2_FS_RS.avi';
% 
% 
% n=n+1;
% files(n).subj = 'Prey13_LT';
% files(n).lighting = OFF;
% files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey_13_LT_trial_3_NIR_LS\Tracks\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_whitebgd_4_15_16\Prey_13_LT_trial_3_NIR_LS.avi';
% 
% 
% n=n+1;
% files(n).subj = 'Prey13_LT';
% files(n).lighting = ON;
% files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey_13_LT_trial_4_FS_LS\Tracks\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_whitebgd_4_15_16\Prey_13_LT_trial_4_FS_LS.avi';
% 
% 
% n=n+1;
% files(n).subj = 'Prey14_NT';
% files(n).lighting = OFF;
% files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey_14_NT_trial_1_NIR_LS\Tracks\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_whitebgd_4_15_16\Prey_14_NT_trial_1_NIR_LS.avi';
% 
% 
% n=n+1;
% files(n).subj = 'Prey14_NT';
% files(n).lighting = ON;
% files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey_14_NT_trial_2_FS_LS\Tracks\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_whitebgd_4_15_16\Prey_14_NT_trial_2_FS_LS.avi';
% 
% 
% n=n+1;
% files(n).subj = 'Prey14_NT';
% files(n).lighting = ON;
% files(n).trackpts = 'Plexiglass_whitebgd_4_15_16\Prey_14_NT_trial_4_FS_RS\Tracks\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_whitebgd_4_15_16\Prey_14_NT_trial_4_FS_RS.avi';
% 
% 
% n=n+1;
% files(n).subj = 'Prey12_PT';
% files(n).lighting = OFF;
% files(n).trackpts = 'Plexiglass_5_31_16\Prey12_PT_trial3_NIR\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_5_31_16\Prey_12_PT_trial_3_RS_NIR.avi';
% 
% 
% n=n+1;
% files(n).subj = 'Prey13_NT';
% files(n).lighting = OFF;
% files(n).trackpts = 'Plexiglass_5_31_16\Prey13_NT_trial4_NIR\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).Moviefile='Plexiglass_5_31_16\Prey_13_NT_trial_4_RS_NIR.avi';
% 
% 
% n=n+1;
% files(n).subj = 'Prey13_PT';
% files(n).lighting = OFF;
% files(n).trackpts = 'Plexiglass_5_31_16\Prey13_PT_trial4_NIR\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_5_31_16\Prey_13_PT_trial_4_LS_NIR.avi';
% 
% 
% n=n+1;
% files(n).subj = 'Prey14_PT';
% files(n).lighting = OFF;
% files(n).trackpts = 'Plexiglass_5_31_16\Prey_14_PT_trial1_NIR\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_5_31_16\Prey_14_PT_trial_1_RS_NIR.avi';
% 
% 
% n=n+1;
% files(n).subj = 'Prey14_PT';
% files(n).lighting = OFF;
% files(n).trackpts = 'Plexiglass_5_31_16\Prey14_PT_trial3_NIR\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_5_31_16\Prey_14_PT_trial_3_LS_NIR.avi';
% 
% 
% n=n+1;
% files(n).subj = 'Prey12_PT';
% files(n).lighting = OFF;
% files(n).trackpts = 'Plexiglass_5_31_16\Prey12_PT_trail1_NIR_2\DLTdv5_data_xypts.csv';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=60;
% files(n).scale=35;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='Plexiglass_5_31_16\Prey_14_PT_trial_3_LS_NIR.avi';
% 


%% 
n=n+1;
files(n).subj = 'Prey12_NT_NIR_RS_trial2'; 
files(n).lighting = OFF;
files(n).trackpts = 'Plexiglass_6_20_16\Prey12_NT_NIR_RS_trial2\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 100; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=1; %assign difference grey or size values to group numbers 
files(n).Moviefile='Plexiglass_6_20_16\Prey12_NT_NIR_RS_trial2.avi';

%
n=n+1;
files(n).subj = 'Prey12_NT_FS_LS_trial3'; 
files(n).lighting = ON;
files(n).trackpts = 'Plexiglass_6_20_16\Prey12_NT_FS_LS_trial3\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 100; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=1; %assign difference grey or size values to group numbers 
files(n).Moviefile='Plexiglass_6_20_16\Prey12_NT_FS_LS_trial3.avi';
%
n=n+1;
files(n).subj = 'Prey12_NT_FS_RS_trial1'; 
files(n).lighting = ON;
files(n).trackpts = 'Plexiglass_6_20_16\Prey12_NT_FS_RS_trial1\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 100; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=1; %assign difference grey or size values to group numbers 
files(n).Moviefile='Plexiglass_6_20_16\Prey12_NT_FS_RS_trial1.avi';

%
n=n+1;
files(n).subj = 'Prey12_NT_NIR_LS_trial4'; 
files(n).lighting = OFF;
files(n).trackpts = 'Plexiglass_6_20_16\Prey12_NT_NIR_LS_Trial4\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 100; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=1; %assign difference grey or size values to group numbers 
files(n).Moviefile='Plexiglass_6_20_16\Prey12_NT_NIR_LS_Trial4.avi';

%
n=n+1;
files(n).subj = 'Prey12_pT_FS_RS_trial2'; 
files(n).lighting = ON;
files(n).trackpts = 'Plexiglass_6_20_16\Prey12_pT_FS_RS_trial2\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 100; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=1; %assign difference grey or size values to group numbers 
files(n).Moviefile='Plexiglass_6_20_16\Prey12_pT_FS_RS_trial2.avi';


%% Emily's size virtual stim data
n=n+1;
files(n).subj = 'Prey20_TT_slide2'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_TT_size_slide2\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 4; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=3; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_TT_size_slide2.avi';


n=n+1;
files(n).subj = 'Prey20_TT_slide3'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey 20_TT_ size_slide3\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 2; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=2; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_TT_size_slide2.avi';


n=n+1;
files(n).subj = 'Prey20_TT_size_slide4'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_TT_size_slide4\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 1; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=1; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_TT_size_slide4.avi';

n=n+1;
files(n).subj = 'Prey20_TT_size_slide5'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_TT_size_slide5\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.5; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=4; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_TT_size_slide5.avi';

n=n+1;
files(n).subj = 'Prey20_TT_size_slide6'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_TT_size_slide6\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.25; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=5; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_TT_size_slide6.avi';

n=n+1;
files(n).subj = 'Prey20_TT_size_slide7'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_TT_size_slide7\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.125; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=6; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_TT_size_slide7.avi';

n=n+1;
files(n).subj = 'Prey20_RT_size_slide9'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_RT_size_slide9\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.5; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=4; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_RT_size_slide9.avi';

n=n+1;
files(n).subj = 'Prey20_RT_size_slide10'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_RT_size_slide10\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 2; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=2; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_RT_size_slide10.avi';

n=n+1;
files(n).subj = 'Prey20_RT_size_slide11'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_RT_size_slide11\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 1; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=1; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_RT_size_slide11.avi';

n=n+1;
files(n).subj = 'Prey20_RT_size_slide 12'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_RT_size_slide 12\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 2; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=2; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_RT_size_slide 12.avi';

n=n+1;
files(n).subj = 'Prey20_RT_size_slide13'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_RT_size_slide13\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.125; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=6; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_RT_size_slide13.avi';

n=n+1;
files(n).subj = 'Prey20_RT_size_slide14'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_RT_size_slide14\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.25; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=5; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_RT_size_slide14.avi';

n=n+1;
files(n).subj = 'Prey20_LT_size_slide16'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_LT_size_slide16\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.5; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=4; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_LT_size_slide16.avi';

n=n+1;
files(n).subj = 'Prey20_LT_size_slide17'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_LT_size_slide17\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 2; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=2; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_LT_size_slide17.avi';

% n=n+1;
% files(n).subj = 'Prey20_LT_size_slide18'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_LT_size_slide18\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 1; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=1; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_LT_size_slide18.avi';

n=n+1;
files(n).subj = 'Prey20_LT_size_slide19'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_LT_size_slide19\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 4; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=3; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_LT_size_slide19.avi';

n=n+1;
files(n).subj = 'Prey20_LT_size_slide20'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_LT_size_slide20\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.125; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=6; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_LT_size_slide20.avi';

n=n+1;
files(n).subj = 'Prey20_LT_size_slide21'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_LT_size_slide21\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.25; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=5; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_LT_size_slide21.avi';

n=n+1;
files(n).subj = 'Prey20_PT2(big)_size_slide23'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_PT2(big)_size_slide23\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.125; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=6; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_PT2(big)_size_slide23.avi';

n=n+1;
files(n).subj = 'Prey20_PT2(big)_size_slide24'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_PT2(big)_size_slide24\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 2; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=2; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_PT2(big)_size_slide24.avi';

n=n+1;
files(n).subj = 'Prey20_PT2(big)_size_slide25'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_PT2(big)_size_slide25\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.5; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=4; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_PT2(big)_size_slide25.avi';

n=n+1;
files(n).subj = 'Prey20_PT2(big)_size_slide26'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_PT2(big)_size_slide26\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.25; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=5; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_PT2(big)_size_slide26.avi';

n=n+1;
files(n).subj = 'Prey20_PT2(big)_size_slide27'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_PT2(big)_size_slide27\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 1; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=1; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_PT2(big)_size_slide27.avi';

n=n+1;
files(n).subj = 'Prey20_PT2(big)_size_slide28'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_PT2(big)_size_slide28\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 4; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=3; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_PT2(big)_size_slide28.avi';

n=n+1;
files(n).subj = 'Prey20_PT1_size_slide30'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_PT1_size_slide30\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.125; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=6; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_PT1_size_slide30.avi';

n=n+1;
files(n).subj = 'Prey20_PT1_size_slide31'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_PT1_size_slide31\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 2; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=2; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_PT1_size_slide31.avi';


n=n+1;
files(n).subj = 'Prey20_PT1_size_slide32_2'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_PT1_size_slide32_2\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.25; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=5; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_PT1_size_slide32_2.avi';

n=n+1;
files(n).subj = 'Prey20_PT1_size_slide33'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_PT1_size_slide33\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.5; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=4; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_PT1_size_slide33.avi';

n=n+1;
files(n).subj = 'Prey20_PT1_size_slide34'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_PT1_size_slide34\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 4; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=3; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_PT1_size_slide34.avi';

n=n+1;
files(n).subj = 'Prey20_PT1_size_slide35'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_virtualStim_8_11_16\Prey20_PT1_size_slide35\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 1; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=1; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_virtualStim_8_11_16\Prey20_PT1_size_slide35.avi';

%% new day's worth of data
n=n+1;
files(n).subj = 'Prey24_PT_slide2'; 
files(n).lighting = ON;
files(n).trackpts = 'Size_9_1_16\Prey24_PT_slide2\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 4; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=3; %assign difference grey or size values to group numbers 
files(n).Moviefile='Size_9_1_16\Prey24_PT_slide2.avi';

% n=n+1;
% files(n).subj = 'Prey24_PT_slide3'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_PT_slide3\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 2; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=2; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_PT_slide3.avi';

% n=n+1;
% files(n).subj = 'Prey24_PT_slide4'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_PT_slide4\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 1; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=1; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_PT_slide4.avi';
% 
% n=n+1;
% files(n).subj = 'Prey24_PT_slide5'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_PT_slide5\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 0.5; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=4; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_PT_slide5.avi';
% 
% n=n+1;
% files(n).subj = 'Prey24_PT_slide6'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_PT_slide6\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 0.25; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=5; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_PT_slide6.avi';
% 
% n=n+1;
% files(n).subj = 'Prey24_PT_slide7'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_PT_slide7\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 0.125; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=6; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_PT_slide7.avi';
% 
% n=n+1;
% files(n).subj = 'Prey24_NT_slide9'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_NT_slide9\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 0.5; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=4; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_NT_slide9.avi';
% 
% n=n+1;
% files(n).subj = 'Prey 24_NT_slide10'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey 24_NT_slide10\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 2; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=2; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey 24_NT_slide10.avi';
% 
% n=n+1;
% files(n).subj = 'Prey24_NT_slide11'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_NT_slide11\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 1; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=1; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_NT_slide11.avi';
% 
% n=n+1;
% files(n).subj = 'Prey24_NT_slide13'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_NT_slide13\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 0.125; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=6; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_NT_slide13.avi';
% 
% n=n+1;
% files(n).subj = 'Prey24_NT_slide14'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_NT_slide14\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 0.25; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=5; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_NT_slide14.avi';
% 
% n=n+1;
% files(n).subj = 'Prey24_TT_slide16'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_TT_slide16\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 0.5; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=4; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_TT_slide16.avi';
% 
% n=n+1;
% files(n).subj = 'Prey24_TT_slide17'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_TT_slide17\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 2; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=2; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_TT_slide17.avi';
% 
% n=n+1;
% files(n).subj = 'Prey24_TT_slide18'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_TT_slide18\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 1; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=1; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_TT_slide18.avi';
% 
% n=n+1;
% files(n).subj = 'Prey24_TT_slide19'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_TT_slide19\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 4; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=3; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_TT_slide19.avi';
% 
% % n=n+1;
% % files(n).subj = 'Prey24_TT_slide20'; 
% % files(n).lighting = ON;
% % files(n).trackpts = 'Size_9_1_16\Prey24_TT_slide20\DLTdv5_data_xypts.csv';%specific path and file name
% % files(n).contrast = 100; %enter grey filter ID
% % files(n).size = 0.125; % enter scaled virtual stimuli ID
% % files(n).notes = '';
% % files(n).fps=60; 
% % files(n).scale=28;%pixels/cm on video tracking
% % files(n).group=6; %assign difference grey or size values to group numbers 
% % files(n).Moviefile='Size_9_1_16\Prey24_TT_slide20.avi';
% % 
% n=n+1;
% files(n).subj = 'Prey24_TT_slide21'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_TT_slide21\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 0.25; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=5; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_TT_slide21.avi';

% n=n+1;
% files(n).subj = 'Prey24_RT_slide23'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_RT_slide23\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 0.125; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=6; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_RT_slide23.avi';
% 
% n=n+1;
% files(n).subj = 'Prey24_RT_slide24'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_RT_slide24\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 2; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=2; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_RT_slide24.avi';
% 
% n=n+1;
% files(n).subj = 'Prey24_RT_slide25'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_RT_slide25\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 0.5; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=4; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_RT_slide25.avi';
% 
% n=n+1;
% files(n).subj = 'Prey24_RT_slide26'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_RT_slide26\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 0.25; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=5; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_RT_slide26.avi';
% 
% n=n+1;
% files(n).subj = 'Prey24_RT_slide27'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_RT_slide27\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 1; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=1; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_RT_slide27.avi';
% 
% n=n+1;
% files(n).subj = 'Prey24_RT_slide28'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_RT_slide28\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 4; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=3; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_RT_slide28.avi';
% 
% n=n+1;
% files(n).subj = 'Prey24_LT_slide30'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_LT_slide30\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 0.125; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=6; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_LT_slide30.avi';
% 
% n=n+1;
% files(n).subj = 'Prey24_LT_slide31'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_LT_slide31\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 2; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=2; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_LT_slide31.avi';
% 
% 
% n=n+1;
% files(n).subj = 'Prey24_LT_slide32'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_LT_slide32\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 0.5; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=4; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_LT_slide32.avi';
%
% n=n+1;
% files(n).subj = 'Prey24_LT_slide33'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_LT_slide33\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 1; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=1; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_LT_slide33.avi';
%
% n=n+1;
% files(n).subj = 'Prey24_LT_slide34'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_LT_slide34\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 4; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=3; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_LT_slide34.avi';
%
% n=n+1;
% files(n).subj = 'Prey24_LT_slide35'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Size_9_1_16\Prey24_LT_slide35\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 1; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=1; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Size_9_1_16\Prey24_LT_slide35.avi';



%% Dolly's size virtual stim, naive

n=n+1;
files(n).subj = 'Prey21_BL-LT_size_slide9'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey21_BL-LT_size_slide9\tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.5; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=10; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey21_BL-LT_size_slide9.avi';
%
n=n+1;
files(n).subj = 'Prey21_BL-LT_size_slide10'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey21_BL-LT_size_slide10\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 4; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=9; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey21_BL-LT_size_slide10.avi';
%
n=n+1;
files(n).subj = 'Prey21_BL-LT_size_slide11'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey21_BL-LT_size_slide11\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 1; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=7; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey21_BL-LT_size_slide11.avi';
%
% n=n+1;
% files(n).subj = 'Prey21_BL-LT_size_slide12'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Naive_mice\8_30_16\Prey21_BL-LT_size_slide12\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 2; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=8; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Naive_mice\8_30_16\Prey21_BL-LT_size_slide12.avi';
%
n=n+1;
files(n).subj = 'Prey21_BL-LT_size_slide13'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey21_BL-LT_size_slide13\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.125; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=12; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey21_BL-LT_size_slide13.avi';
%
n=n+1;
files(n).subj = 'Prey21_BL-LT_size_slide14'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey21_BL-LT_size_slide14\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.25; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=11; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey21_BL-LT_size_slide14.avi';
%
n=n+1;
files(n).subj = 'Prey21_BR-LT_size_slide2'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey21_BR-LT_size_slide2\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 4; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=9; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey21_BR-LT_size_slide2.avi';
%
n=n+1;
files(n).subj = 'Prey21_BR-LT_size_slide3'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey21_BR-LT_size_slide3\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 2; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=8; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey21_BR-LT_size_slide3.avi';
%
n=n+1;
files(n).subj = 'Prey21_BR-LT_size_slide4'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey21_BR-LT_size_slide4\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 1; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=7; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey21_BR-LT_size_slide4.avi';
%
n=n+1;
files(n).subj = 'Prey21_BR-LT_size_slide5'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey21_BR-LT_size_slide5\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.5; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=10; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey21_BR-LT_size_slide5.avi';
%
n=n+1;
files(n).subj = 'Prey21_BR-LT_size_slide6'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey21_BR-LT_size_slide6\tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.25; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=11; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey21_BR-LT_size_slide6.avi';
%
n=n+1;
files(n).subj = 'Prey21_BR-LT_size_slide7'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey21_BR-LT_size_slide7\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.125; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=12; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey21_BR-LT_size_slide7.avi';
%
n=n+1;
files(n).subj = 'Prey22_NT_Nt_size_slide16'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey22_NT_Nt_size_slide16\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.5; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=10; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey22_NT_Nt_size_slide16.avi';
%
n=n+1;
files(n).subj = 'Prey22_NT_Nt_size_slide17'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey22_NT_Nt_size_slide17\TRACKS\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 2; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=8; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey22_NT_Nt_size_slide17.avi';
%
n=n+1;
files(n).subj = 'Prey22_NT_Nt_size_slide18'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey22_NT_Nt_size_slide18\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 1; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=7; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey22_NT_Nt_size_slide18.avi';
%
n=n+1;
files(n).subj = 'Prey22_NT_Nt_size_slide19'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey22_NT_Nt_size_slide19\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 4; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=9; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey22_NT_Nt_size_slide19.avi';
%
n=n+1;
files(n).subj = 'Prey22_NT_Nt_size_slide20'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey22_NT_Nt_size_slide20\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.125; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=12; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey22_NT_Nt_size_slide20.avi';
%
n=n+1;
files(n).subj = 'Prey22_NT_Nt_size_slide21'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey22_NT_Nt_size_slide21\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.25; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=11; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey22_NT_Nt_size_slide21.avi';
%
n=n+1;
files(n).subj = 'Prey22_NT-pstripe_size_slide23'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey22_NT-pstripe_size_slide23\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.125; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=12; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey22_NT-pstripe_size_slide23.avi';
%
n=n+1;
files(n).subj = 'Prey22_NT-pstripe_size_slide24'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey22_NT-pstripe_size_slide24\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 2; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=8; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey22_NT-pstripe_size_slide24.avi';
%
n=n+1;
files(n).subj = 'Prey22_NT-pstripe_size_slide25'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey22_NT-pstripe_size_slide25\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.5; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=10; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey22_NT-pstripe_size_slide25.avi';
%
n=n+1;
files(n).subj = 'Prey22_NT-pstripe_size_slide26'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey22_NT-pstripe_size_slide26\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.25; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=11; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey22_NT-pstripe_size_slide26.avi';
%
n=n+1;
files(n).subj = 'Prey22_NT-pstripe_size_slide27'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey22_NT-pstripe_size_slide27\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 1; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=7; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey22_NT-pstripe_size_slide27.avi';
%
n=n+1;
files(n).subj = 'Prey22_NT-pstripe_size_slide28'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey22_NT-pstripe_size_slide28\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 4; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=9; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey22_NT-pstripe_size_slide28.avi';
%
n=n+1;
files(n).subj = 'Prey22_NT-pt_size_slide30'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey22_NT-pt_size_slide30\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.125; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=12; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey22_NT-pt_size_slide30.avi';
%
n=n+1;
files(n).subj = 'Prey22_NT-pt_size_slide31'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey22_NT-pt_size_slide31\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 2; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=8; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey22_NT-pt_size_slide31.avi';
%
n=n+1;
files(n).subj = 'Prey22_NT-pt_size_slide32'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey22_NT-pt_size_slide32\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.25; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=11; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey22_NT-pt_size_slide32.avi';
%
n=n+1;
files(n).subj = 'Prey22_NT-pt_size_slide33'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey22_NT-pt_size_slide33\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.5; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=10; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey22_NT-pt_size_slide33.avi';
%
n=n+1;
files(n).subj = 'Prey22_NT-pt_size_slide34'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey22_NT-pt_size_slide34\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 4; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=9; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey22_NT-pt_size_slide34.avi';
%
n=n+1;
files(n).subj = 'Prey22_NT-pt_size_slide354'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\8_30_16\Prey22_NT-pt_size_slide354\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 1; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=7; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\8_30_16\Prey22_NT-pt_size_slide354.avi';
%
n=n+1;
files(n).subj = 'Pre21_BK-RT_size_slide9'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Pre21_BK-RT_size_slide9\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.5; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=10; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Pre21_BK-RT_size_slide9.avi';
%
n=n+1;
files(n).subj = 'Pre21_BK-RT_size_slide10'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Pre21_BK-RT_size_slide10\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 4; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=9; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Pre21_BK-RT_size_slide10.avi';
%
n=n+1;
files(n).subj = 'Pre21_BK-RT_size_slide11'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Pre21_BK-RT_size_slide11\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 1; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=7; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Pre21_BK-RT_size_slide11.avi';
%
n=n+1;
files(n).subj = 'Pre21_BK-RT_size_slide12'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Pre21_BK-RT_size_slide12\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 2; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=8; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Pre21_BK-RT_size_slide12.avi';
%
% n=n+1;
% files(n).subj = 'Pre21_BK-RT_size_slide13'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Naive_mice\9_13_16\Pre21_BK-RT_size_slide13\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; %enter grey filter ID
% files(n).size = 0.125; % enter scaled virtual stimuli ID
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=12; %assign difference grey or size values to group numbers 
% files(n).Moviefile='Naive_mice\9_13_16\Pre21_BK-RT_size_slide13.avi';
% %
n=n+1;
files(n).subj = 'Pre21_BK-RT_size_slide14'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Pre21_BK-RT_size_slide14\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.25; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=11; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Pre21_BK-RT_size_slide14.avi';
%
n=n+1;
files(n).subj = 'Prey21_BK-TT_size_slide23'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Prey21_BK-TT_size_slide23\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.125; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=12; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Prey21_BK-TT_size_slide23.avi';
%
n=n+1;
files(n).subj = 'Prey21_BK-TT_size_slide24'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Prey21_BK-TT_size_slide24\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 2; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=8; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Prey21_BK-TT_size_slide24.avi';
%
n=n+1;
files(n).subj = 'Prey21_BK-TT_size_slide25'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Prey21_BK-TT_size_slide25\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.5; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=10; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Prey21_BK-TT_size_slide25.avi';
%
n=n+1;
files(n).subj = 'Prey21_BK-TT_size_slide26'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Prey21_BK-TT_size_slide26\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.25; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=11; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Prey21_BK-TT_size_slide26.avi';
%
n=n+1;
files(n).subj = 'Prey21_BK-TT_size_slide27'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Prey21_BK-TT_size_slide27\tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 1; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=7; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Prey21_BK-TT_size_slide27.avi';
%
n=n+1;
files(n).subj = 'Prey21_BK-TT_size_slide28'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Prey21_BK-TT_size_slide28\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 4; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=9; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Prey21_BK-TT_size_slide28.avi';
%
n=n+1;
files(n).subj = 'Prey21_BR-RT_size_slide16'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Prey21_BR-RT_size_slide16\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.5; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=10; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Prey21_BR-RT_size_slide16.avi';
%
n=n+1;
files(n).subj = 'Prey21_BR-RT_size_slide17'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Prey21_BR-RT_size_slide17\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 2; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=8; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Prey21_BR-RT_size_slide17.avi';
%
n=n+1;
files(n).subj = 'Prey21_BR-RT_size_slide18'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Prey21_BR-RT_size_slide18\tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 1; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=7; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Prey21_BR-RT_size_slide18.avi';
%
n=n+1;
files(n).subj = 'Prey21_BR-RT_size_slide19'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Prey21_BR-RT_size_slide19\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 4; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=9; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Prey21_BR-RT_size_slide19.avi';
%
n=n+1;
files(n).subj = 'Prey21_BR-RT_size_slide20'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Prey21_BR-RT_size_slide20\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.125; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=12; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Prey21_BR-RT_size_slide20.avi';
%
n=n+1;
files(n).subj = 'Prey21_BR-RT_size_slide21'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Prey21_BR-RT_size_slide21\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.25; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=11; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Prey21_BR-RT_size_slide21.avi';
%
n=n+1;
files(n).subj = 'Prey22_RT_size_slide3'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Prey22_RT_size_slide3\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 2; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=8; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Prey22_RT_size_slide3.avi';
%
n=n+1;
files(n).subj = 'Prey22_RT_size_slide4'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Prey22_RT_size_slide4\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 1; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=7; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Prey22_RT_size_slide4.avi';
%
n=n+1;
files(n).subj = 'Prey22_RT_size_slide5'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Prey22_RT_size_slide5\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.5; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=10; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Prey22_RT_size_slide5.avi';
%
n=n+1;
files(n).subj = 'Prey22_RT_size_slide6'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Prey22_RT_size_slide6\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.25; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=11; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Prey22_RT_size_slide6.avi';
%
n=n+1;
files(n).subj = 'Prey22_RT_size_slide7'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Prey22_RT_size_slide7\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 0.125; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=12; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Prey22_RT_size_slide7.avi';
%
n=n+1;
files(n).subj = 'Prey22_RT_size_slide2'; 
files(n).lighting = ON;
files(n).trackpts = 'Naive_mice\9_13_16\Prey22_RT_size_slide2\Tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; %enter grey filter ID
files(n).size = 4; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=9; %assign difference grey or size values to group numbers 
files(n).Moviefile='Naive_mice\9_13_16\Prey22_RT_size_slide2.avi';
%
