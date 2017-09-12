ON=1; OFF=0;
pathname = 'X:\Preycapture\';%miyazaki path to prey capture folder stored on Eyre 
n=0;

%%%group:1=pre-ctl,2=pre-FM, 3=post-ctl, 4=post-FM

%%%female=1 and male=2;
%%%contrast:100
%%%size:100
%%
n=n+1;
files(n).subj = '9_29_16_mouseA';
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\9_29_16\MouseA_virtual_1\DLTdv5_data_xypts.csv';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=28;%pixels/cm on video tracking
files(n).group=4;%muscimol condition
files(n).Moviefile='';

% 
n=n+1;
files(n).subj = '10_25_16_MouseA_pinkstrip_virtual_trial1'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\10_25_16\MouseA_pinkstrip_virtual_trial1\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=3; %control condition
files(n).Moviefile='';
% 
n=n+1;
files(n).subj = '10_25_16_mouseA_virtual_trial2'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\10_25_16\MouseA_pinkstrip_virtual_trial2\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=3; %condition of injections
files(n).Moviefile='';
% % 
n=n+1;
files(n).subj = '10_25_16_MouseB_pinktip_virtual_trial1'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\10_25_16\MouseB_pinktip_virtual_trial1\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=4; %condition of injections
files(n).Moviefile='';
% % 
n=n+1;
files(n).subj = '10_25_16_MouseB_pinktip_virtual_trial2'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\10_25_16\MouseB_pinktip_virtual_trial2\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=4; %condition of injections
files(n).Moviefile='';
% % 
n=n+1;
files(n).subj = '10_25_16_pinkstrip_pre_virtual'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\10_25_16\pinkstrip_pre_virtual\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=1; %condition of injections
files(n).Moviefile='';
% 
% n=n+1;
% files(n).subj = '10_25_16_pinkstrip_pre_virtual'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Muscimol\10_25_16\pinkstrip_pre_virtual\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; 
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=1; %condition of injections
% files(n).Moviefile='';
% 
n=n+1;
files(n).subj = '10_25_16_pinktip_pre_virtual'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\10_25_16\Pinktip_pre_virtual\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=2; %condition of injections
files(n).Moviefile='';
% 
n=n+1;
files(n).subj = '10_27_16_mouseA_PT_large_virtual_trial1'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\10_27_16\mouseA_PT_large_virtual_trial1\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=4; %condition of injections
files(n).Moviefile='';
% 
% n=n+1;
% files(n).subj = '10_27_16_mouseA_PT_large_virtual_trial2'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Muscimol\10_27_16\mouseA_PT_large_virtual_trial2\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; 
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=4; %condition of injections
% files(n).Moviefile='';
% not tracked correctly

n=n+1;
files(n).subj = '10_27_16_mouseB_blacktip_virtual_trial1'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\10_27_16\mouseB_blacktip_virtual_trial1\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=3; %condition of injections
files(n).Moviefile='';
% % 
n=n+1;
files(n).subj = '11_1_16_MouseA_virtual_pre_musc'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_1_16\MouseA_virtual_pre_musc\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=1; %condition of injections
files(n).Moviefile='';
% 
n=n+1;
files(n).subj = '11_1_16_MouseA_virtual_trial1'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_1_16\MouseA_virtual_trial1\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=3; %condition of injections
files(n).Moviefile='';
% 
n=n+1;
files(n).subj = '11_1_16_MouseA_virtual_trial2'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_1_16\MouseA_virtual_trial2\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=3; %condition of injections
files(n).Moviefile='';
% % 
n=n+1;
files(n).subj = '11_1_16_MouseB_virtual_pre_musc'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_1_16\MouseB_virtual_pre_musc\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=2; %condition of injections
files(n).Moviefile='';
% % 
n=n+1;
files(n).subj = '11_1_16_MouseB_virtual_trial1'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_1_16\MouseB_virtual_trial1\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=4; %condition of injections
files(n).Moviefile='';
% % 
n=n+1;
files(n).subj = '11_1_16_MouseB_virtual_trial2'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_1_16\MouseB_virtual_trial2\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=4; %condition of injections
files(n).Moviefile='';
% % 
n=n+1;
files(n).subj = '11_1_16_MouseC_virtual_pre_musc'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_1_16\MouseC_virtual_pre_musc\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=2; %condition of injections
files(n).Moviefile='';
% 
% n=n+1;
% files(n).subj = '11_1_16_MouseC_virtual_trial1'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Muscimol\11_1_16\MouseC_virtual_trial1\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; 
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=4; %condition of injections
% files(n).Moviefile='';
% no data in excel file, not tracked?
n=n+1;
files(n).subj = '11_17_16_MouseA_pre_virtual'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_17_16\MouseA_pre_virtual\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=2; %condition of injections
files(n).Moviefile='';
% 
n=n+1;
files(n).subj = '11_17_16_MouseA_virtual_trial1'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_17_16\MouseA_virtual_trial1\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=4; %condition of injections
files(n).Moviefile='';
% 
n=n+1;
files(n).subj = '11_17_16_MouseA_virtual_trial2'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_17_16\MouseA_virtual_trial2\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=4; %condition of injections
files(n).Moviefile='';
% 
n=n+1;
files(n).subj = '11_17_16_MouseB_virtual_trial1'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_17_16\MouseB_virtual_trial1\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=4; %condition of injections
files(n).Moviefile='';
% 
n=n+1;
files(n).subj = '11_17_16_MouseB_virtual_trial2'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_17_16\MouseB_virtual_trial2\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=4; %condition of injections
files(n).Moviefile='';
% 
n=n+1;
files(n).subj = '11_17_16_MouseB3_pre_virtual'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_17_16\MouseB3_pre_virtual\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=2; %condition of injections
files(n).Moviefile='';
% 
n=n+1;
files(n).subj = '11_17_16_MouseC_pre_virtual'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_17_16\MouseC_pre_virtual\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=1; %condition of injections
files(n).Moviefile='';
%
n=n+1;
files(n).subj = '11_17_16_MouseC_virtual_trial1'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_17_16\MouseC_virtual_trial1\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=3; %condition of injections
files(n).Moviefile='';
%
n=n+1;
files(n).subj = '11_17_16_MouseC_virtual_trial2'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_17_16\MouseC_virtual_trial2\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=3; %condition of injections
files(n).Moviefile='';
%
%
%% %%%group:1=pre-ctl,2=pre-FM, 3=post-ctl, 4=post-FM
n=n+1;
files(n).subj = '11_22_16_mouse1_pre_virtual'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_22_16\mouse1_pre_virtual_tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=2; %condition of injections
files(n).Moviefile='';
%
n=n+1;
files(n).subj = '11_22_16_mouse1_virtual_trial1'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_22_16\mouse1_virtual_trial1\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=4; %condition of injections
files(n).Moviefile='';
%
n=n+1;
files(n).subj = '11_22_16_mouse1_virtual_trial2'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_22_16\mouse1_virtual_trial2\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=4; %condition of injections
files(n).Moviefile='';
%
% n=n+1;
% files(n).subj = '11_22_16_mouse2_pre_virtual'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Muscimol\11_22_16\mouse2_pre_virtual_tracks\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; 
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=0; %condition of injections
% files(n).Moviefile='';
% %
% n=n+1;
% files(n).subj = '11_22_16_mouse2_virtual_trial1'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Muscimol\11_22_16\mouse2_virtual_trial1_tracks\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; 
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=0; %condition of injections
% files(n).Moviefile='';
% %
% n=n+1;
% files(n).subj = '11_22_16_mouse2_virtual_trial2'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Muscimol\11_22_16\mouse2_virtual_trial2_tracks\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; 
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=0; %condition of injections
% files(n).Moviefile='';
% %%%%group:1=pre-ctl,2=pre-FM, 3=post-ctl, 4=post-F
n=n+1;
files(n).subj = '11_22_16_mouse3_pre_virtual'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_22_16\mouse3_pre_virtual_tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=2; %condition of injections
files(n).Moviefile='';
%
n=n+1;
files(n).subj = '11_22_16_mouse3_virtual_trial1'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_22_16\mouse3_virtual_trial1_tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=4; %condition of injections
files(n).Moviefile='';
%
n=n+1;
files(n).subj = '11_22_16_mouse3_virtual_trial2'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_22_16\mouse3_virtual_trial2_tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=4; %condition of injections
files(n).Moviefile='';
%%
n=n+1;
files(n).subj = '11_22_16_mouse4_pre_virtual'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_22_16\mouse4_pre_virtual_tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=1; %condition of injections
files(n).Moviefile='';
%
n=n+1;
files(n).subj = '11_22_16_mouse4_virtual_trial1'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\11_22_16\mouse4_virtual_trial1_tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=3; %condition of injections
files(n).Moviefile='';
%
 
%% enter DREADD data
% group1=negative virus; group2= Postive for virus 

%
n=n+1;
files(n).subj = '4_12_17\NTSR1_1_3_LN_virtual_1X_stim_left_side_'; 
files(n).lighting = ON;
files(n).trackpts = 'DREADDS_Tracks\LN_virtual_1x_stim_left_side_tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=30; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=5; %condition of injections
files(n).Moviefile='';
%
n=n+1;
files(n).subj = '4_12_17\NTSR1_1_3_LN\virtual_pt5X_stim_left_side'; 
files(n).lighting = ON;
files(n).trackpts = 'DREADDS_Tracks\LN_virtual_pt5x_stim_left_side\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=30; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=5; %condition of injections
files(n).Moviefile='';
%
n=n+1;
files(n).subj = '4_12_17\NTSR1_1_3_LN\virtual_pt5X_stim_middle'; 
files(n).lighting = ON;
files(n).trackpts = 'DREADDS_Tracks\LN_virtual_pt5x_stim_middle_tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=30; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=5; %condition of injections
files(n).Moviefile='';
%
n=n+1;
files(n).subj = '4_12_17\NTSR1_1_3_LN\virtual_pt5X_stim_right_side'; 
files(n).lighting = ON;
files(n).trackpts = 'DREADDS_Tracks\LN_virtual_pt5x_stim_right_side_tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=30; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=5; %condition of injections
files(n).Moviefile='';
%
n=n+1;
files(n).subj = '4_12_17\NTSR1_1_4_RT\virtual_middle_pt5x'; 
files(n).lighting = ON;
files(n).trackpts = 'DREADDS_Tracks\RT_virtual_middle_pt5x_tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=30; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=6; %condition of injections
files(n).Moviefile='';
%
n=n+1;
files(n).subj = '4_12_17\NTSR1_1_4_RT\virtual_pt5x_stim_rightside(1)'; 
files(n).lighting = ON;
files(n).trackpts = 'DREADDS_Tracks\RT_virtual_pt5x_stim_rightside(1)_tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=30; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=6; %condition of injections
files(n).Moviefile='';
%
n=n+1;
files(n).subj = '4_12_17\NTSR1_1_4_RT\virtual_pt5x_stim_rightside'; 
files(n).lighting = ON;
files(n).trackpts = 'DREADDS_Tracks\RT_virtual_pt5x_stim_rightside_tracks\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=30; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=6; %condition of injections
files(n).Moviefile='';
%


%% V1 cortex muscimol lessions 


n=n+1;
files(n).subj = '1_13_17_mouseA'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol V1 Saved Tracks\1_13_17\MouseA_virtual_trial1\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=7; %condition of injections
files(n).Moviefile='';

n=n+1;
files(n).subj = '1_13_17_mouseA'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol V1 Saved Tracks\1_13_17\MouseA_virtual_trial2\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=7; %condition of injections
files(n).Moviefile='';

n=n+1;
files(n).subj = '1_18_17_mouseB'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol V1 Saved Tracks\1_18_17\MouseB_post_virtual_trial1\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=7; %condition of injections
files(n).Moviefile='';

n=n+1;
files(n).subj = '1_18_17_mouseB'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol V1 Saved Tracks\1_18_17\MouseB_post_virtual_trial2\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=7; %condition of injections
files(n).Moviefile='';

n=n+1;
files(n).subj = '1_18_17_mouseB'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol V1 Saved Tracks\1_18_17\MouseB_post_virtual_trial3\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=7; %condition of injections
files(n).Moviefile='';

n=n+1;
files(n).subj = '1_13_17_mouseB'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\V1\1_13_17\MouseB_virtual_trial1\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=7; %condition of injections
files(n).Moviefile='';

% n=n+1;
% files(n).subj = '1_24_17_mouseA'; 
% files(n).lighting = ON;
% files(n).trackpts = 'Muscimol\V1\1_24_17\MouseA_virtual_trial1\DLTdv5_data_xypts.csv';%specific path and file name
% files(n).contrast = 100; 
% files(n).notes = '';
% files(n).fps=60; 
% files(n).scale=28;%pixels/cm on video tracking
% files(n).group=7; %condition of injections
% files(n).Moviefile='';

n=n+1;
files(n).subj = '1_24_17_mouseA'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\V1\1_24_17\MouseA_virtual_trial2\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=7; %condition of injections
files(n).Moviefile='';

n=n+1;
files(n).subj = '2_9_17_mouseA'; 
files(n).lighting = ON;
files(n).trackpts = 'Muscimol\V1\2_9_17\MouseA_virtual_trial2\DLTdv5_data_xypts.csv';%specific path and file name
files(n).contrast = 100; 
files(n).notes = '';
files(n).fps=60; 
files(n).scale=28;%pixels/cm on video tracking
files(n).group=7; %condition of injections
files(n).Moviefile='';
