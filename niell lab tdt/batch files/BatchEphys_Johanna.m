% file entry info for batch processing of single unit analysis E-phys
% data with or without optogenetic tagging

pathname = 'D:\';%computer drivew where data folders exist

n=0;

%%%% 

%n=n+1;
% files(n).condition=1;%control =1;mutant/drug/etc.=2; 
% files(n).sex=2; %1=female and 2=male
% files(n).age=1;%1=adult, 2=p21, 3= EO3, 4=EO1 etc.
% files(n).expt = '1_12_16';%could enter date in order to organize by date
% files(n).notes = 'post_pinp'; %could comment whether experiment was good bad or other flags
% 
% %file parameters and paths
% 
% files(n).path = 'D:\johanna_analysis\1_12_16\';
% files(n).clusterfile = 'cluster_data_1_12_16_cluster.mat';
% files(n).analysisfile = 'analysis.mat';
% files(n).pathtank='D:\Johanna_tanks\1_12_16\';
% 
% %block info
% files(n).blockDrift = {'drift_02'};
% files(n).blockWn = {'wn_02'};
% files(n).blockPinp = {'laser_4_nw_thresh'};
% files(n).blockBar = {''};
% files(n).blockNoStim={''};
% files(n).postBlocks={''};
% files(n).nchan=64; %number of sites
% 
% %parameters for layer analysis
% files(n).tip_loc_1=575;
% files(n).tip_loc_2=550;
% files(n).angle=65; 


n=n+1;
files(n).condition=1;%control =1;mutant/drug/etc.=2; 
files(n).sex=2; %1=female and 2=male
files(n).age=1;%1=adult, 2=p21, 3= EO3, 4=EO1 etc.
files(n).expt = '1_14_16';%could enter date in order to organize by date
files(n).notes = ''; %could comment whether experiment was good bad or other flags

%file parameters and paths

files(n).path = 'D:\johanna_analysis\1_14_16\';
files(n).clusterfile = 'cluster_data_1_14_16_cluster.mat';
files(n).analysisfile = 'analysis.mat';
files(n).pathtank='D:\Johanna_tanks\1_14_16\';

%block info
files(n).blockDrift = {'drift_01'};
files(n).blockWn = {'wn_01'};
files(n).blockPinp = {'pinp_p7_01' 'pinp_P8_02'};
files(n).prefPinp = {'pinp_P8_02'};
files(n).blockBar = {'bars_01'};
files(n).blockNoStim={''};
files(n).postBlocks={''};
files(n).nchan=64; %number of sites

%parameters for layer analysis
files(n).tip_loc_1=525;
files(n).tip_loc_2=500;
files(n).angle=55; 


n=n+1;
files(n).condition=1;%control =1;mutant/drug/etc.=2; 
files(n).sex=2; %1=female and 2=male
files(n).age=1;%1=adult, 2=p21, 3= EO3, 4=EO1 etc.
files(n).expt = '1_19_16';%could enter date in order to organize by date
files(n).notes = ''; %could comment whether experiment was good bad or other flags

%file parameters and paths

files(n).path = 'D:\johanna_analysis\1_19_16\';
files(n).clusterfile = 'cluster_data_1_19_16_cluster.mat';
files(n).analysisfile = 'analysis.mat';
files(n).pathtank='D:\Johanna_tanks\1_19_16\';

%block info
files(n).blockDrift = {'drift_01'};
files(n).blockWn = {'wn_01'};
files(n).blockPinp = {'laser_4_5' 'laser_5_5' 'Laser_7' 'laser_8'};
files(n).prefPinp = {'laser_8'};
files(n).blockBar = {'bars_01'};
files(n).blockNoStim={''};
files(n).postBlocks={''};
files(n).nchan=64; %number of sites

%parameters for layer analysis
files(n).tip_loc_1=525;
files(n).tip_loc_2=500;
files(n).angle=55; 


n=n+1;
files(n).condition=1;%control =1;mutant/drug/etc.=2; 
files(n).sex=1; %1=female and 2=male
files(n).age=1;%1=adult, 2=p21, 3= EO3, 4=EO1 etc.
files(n).expt = '1_21_16';%could enter date in order to organize by date
files(n).notes = ''; %could comment whether experiment was good bad or other flags

%file parameters and paths

files(n).path = [pathname 'johanna_analysis\1_21_16\'];
files(n).clusterfile = 'cluster_data_1_211_16.mat';
files(n).analysisfile = 'analysis.mat';
files(n).pathtank='D:\Johanna_tanks\1_211_16\';

%block info
files(n).blockDrift = {'drift_01'};
files(n).blockWn = {'wn_01'};
files(n).blockPinp = {'laser_P4_01' 'Laser_4_5_02' 'Laser_p5_03' 'Laser_P6_04' 'Laser_P7_06'};
files(n).prefPinp = {'Laser_4_5_02'};
files(n).blockBar = {'bars_011' 'bars_postpinp'};
files(n).blockNoStim={''};
files(n).postBlocks={''};
files(n).nchan=64; %number of sites

%parameters for layer analysis
files(n).tip_loc_1=450;
files(n).tip_loc_2=425;
files(n).angle=50; 


n=n+1;
files(n).condition=1;%control =1;mutant/drug/etc.=2; 
files(n).sex=1; %1=female and 2=male
files(n).age=1;%1=adult, 2=p21, 3= EO3, 4=EO1 etc.
files(n).expt = '1_26_16';%could enter date in order to organize by date
files(n).notes = ''; %could comment whether experiment was good bad or other flags

%file parameters and paths

files(n).path = [pathname 'johanna_analysis\1_26_16\'];
files(n).clusterfile = 'cluster_data_1_26_16_cluster.mat';
files(n).analysisfile = 'analysis.mat';
files(n).pathtank='D:\Johanna_tanks\1_26_16\';

%block info
files(n).blockDrift = {'drift_01'};
files(n).blockWn = {'wn_01'};
files(n).blockPinp = {'Laser_p3_01_real' 'Laser_p4_02' 'Laser_p5_03' 'laser_p6_04'};
files(n).prefPinp = {'Laser_p4_02'};
files(n).blockBar = {'bars_01' 'bars_02'};
files(n).blockNoStim={''};
files(n).postBlocks={''};
files(n).nchan=64; %number of sites

%parameters for layer analysis
files(n).tip_loc_1=450;
files(n).tip_loc_2=425;
files(n).angle=50; 


n=n+1;
files(n).condition=1;%control =1;mutant/drug/etc.=2; 
files(n).sex=2; %1=female and 2=male
files(n).age=1;%1=adult, 2=p21, 3= EO3, 4=EO1 etc.
files(n).expt = '1_28_16';%could enter date in order to organize by date
files(n).notes = ''; %could comment whether experiment was good bad or other flags

%file parameters and paths

files(n).path = [pathname 'johanna_analysis\1_28_16\'];
files(n).clusterfile = 'cluster_data_1_28_16_JT_cluster.mat';
files(n).analysisfile = 'analysis.mat';
files(n).pathtank='D:\Johanna_tanks\1_28_16_JT\';

%block info
files(n).blockDrift = {'drift_01'};
files(n).blockWn = {'wn_01'};
files(n).blockPinp = {'Laser_p3_01' 'laser_p4_02' 'laser_p4_5_03' 'laser_p5_04' 'laser_p5_5_05'};
files(n).prefPinp = {'laser_p5_5_05'};
files(n).blockBar = {'bars_01' 'bars_02'};
files(n).blockNoStim={''};
files(n).postBlocks={''};
files(n).nchan=64; %number of sites

%parameters for layer analysis
files(n).tip_loc_1=450;
files(n).tip_loc_2=425;
files(n).angle=45; 

n=n+1;
files(n).condition=1;%control =1;mutant/drug/etc.=2; 
files(n).sex=1; %1=female and 2=male
files(n).age=1;%1=adult, 2=p21, 3= EO3, 4=EO1 etc.
files(n).expt = '2_2_16';%could enter date in order to organize by date
files(n).notes = ''; %could comment whether experiment was good bad or other flags

%file parameters and paths

files(n).path = [pathname 'johanna_analysis\2_2_16\'];
files(n).clusterfile = 'cluster_data_2_2_16_cluster.mat';
files(n).analysisfile = 'analysis.mat';
files(n).pathtank='D:\Johanna_tanks\2_2_16\';

%block info
files(n).blockDrift = {'drift_01'};
files(n).blockWn = {'wn_01_real'};
files(n).blockPinp = {'laser_p3_01' 'laser_p3_5_02' 'laser_p4_03' 'Laser_p_4_5_04'};
files(n).prefPinp = {'laser_p3_01'};
files(n).blockBar = {'bars_01' 'bars_02'};
files(n).blockNoStim={''};
files(n).postBlocks={''};
files(n).nchan=64; %number of sites

%parameters for layer analysis
files(n).tip_loc_1=475;
files(n).tip_loc_2=450;
files(n).angle=48; 


n=n+1;
files(n).condition=1;%control =1;mutant/drug/etc.=2; 
files(n).sex=1; %1=female and 2=male
files(n).age=1;%1=adult, 2=p21, 3= EO3, 4=EO1 etc.
files(n).expt = '2_9_16';%could enter date in order to organize by date
files(n).notes = ''; %could comment whether experiment was good bad or other flags

%file parameters and paths

files(n).path = [pathname 'johanna_analysis\2_9_16\'];
files(n).clusterfile = 'cluster_data_2_9_16_cluster.mat';
files(n).analysisfile = 'analysis.mat';
files(n).pathtank='D:\Johanna_tanks\2_9_16\';

%block info
files(n).blockDrift = {'drift_01'};
files(n).blockWn = {'wn_01_real'};
files(n).blockPinp = {'Pinp_p2_5_01' 'pinp_p3_02' 'PinP_p3_5_03'};
files(n).prefPinp = {'Pinp_p2_5_01'};
files(n).blockBar = {'bars_01' 'bars_02'};
files(n).blockNoStim={''};
files(n).postBlocks={''};
files(n).nchan=64; %number of sites

%parameters for layer analysis
files(n).tip_loc_1=450;
files(n).tip_loc_2=425;
files(n).angle=45; 

n=n+1;
files(n).condition=1;%control =1;mutant/drug/etc.=2; 
files(n).sex=1; %1=female and 2=male
files(n).age=1;%1=adult, 2=p21, 3= EO3, 4=EO1 etc.
files(n).expt = '2_18_16';%could enter date in order to organize by date
files(n).notes = ''; %could comment whether experiment was good bad or other flags

%file parameters and paths

files(n).path = [pathname 'johanna_analysis\2_18_16\'];
files(n).clusterfile = 'cluster_data_2_18_16_cluster.mat';
files(n).analysisfile = 'analysis.mat';
files(n).pathtank='D:\Johanna_tanks\2_18_16\';

%block info
files(n).blockDrift = {'drift_new_thresh'};
files(n).blockWn = {'wn_new_thresh'};
files(n).blockPinp = {'laser_p3_01' 'laser_p3_5_02' 'laser_p5_03' 'laser_p5_04' 'laser_p5_5_05' };
files(n).prefPinp = {'laser_p3_01'};
files(n).blockBar = {'bars_new_thresh_real'};
files(n).blockNoStim={''};
files(n).postBlocks={''};
files(n).nchan=64; %number of sites

%parameters for layer analysis
files(n).tip_loc_1=450;
files(n).tip_loc_2=425;
files(n).angle=45; 

n=n+1;
files(n).condition=1;%control =1;mutant/drug/etc.=2; 
files(n).sex=1; %1=female and 2=male
files(n).age=1;%1=adult, 2=p21, 3= EO3, 4=EO1 etc.
files(n).expt = '2_23_16';%could enter date in order to organize by date
files(n).notes = ''; %could comment whether experiment was good bad or other flags

%file parameters and paths

files(n).path = [pathname 'johanna_analysis\2_23_16\'];
files(n).clusterfile = 'cluster_data_2_23_16_cluster.mat';
files(n).analysisfile = 'analysis.mat';
files(n).pathtank='D:\Johanna_tanks\2_23_16\';

%block info
files(n).blockDrift = {'drift_01'};
files(n).blockWn = {'wn_01'};
files(n).blockPinp = {'laser_p3_5_01' 'laser_p4_02' 'laser_p4_5_03'};
files(n).prefPinp = {'laser_p4_5_03'};
files(n).blockBar = {'bars_01' 'bars_02'};
files(n).blockNoStim={''};
files(n).postBlocks={''};
files(n).nchan=64; %number of sites

%parameters for layer analysis
files(n).tip_loc_1=450;
files(n).tip_loc_2=425;
files(n).angle=45; 

n=n+1;
files(n).condition=1;%control =1;mutant/drug/etc.=2; 
files(n).sex=1; %1=female and 2=male
files(n).age=1;%1=adult, 2=p21, 3= EO3, 4=EO1 etc.
files(n).expt = '3_4_16';%could enter date in order to organize by date
files(n).notes = ''; %could comment whether experiment was good bad or other flags

%file parameters and paths

files(n).path = [pathname 'johanna_analysis\3_4_16\'];
files(n).clusterfile = 'cluster_data_3_4_16_cluster.mat';
files(n).analysisfile = 'analysis.mat';
files(n).pathtank='D:\Johanna_tanks\3_4_16\';

%block info
files(n).blockDrift = {'drift_01'};
files(n).blockWn = {'wn_01'};
files(n).blockPinp = {'laser_p3_01' 'laser_p2_5_02' 'laser_p3_5_03' 'laser_p4_4' };
files(n).prefPinp = {'laser_p3_01'};
files(n).blockBar = {'bars_01' 'bars_02'};
files(n).blockNoStim={''};
files(n).postBlocks={''};
files(n).nchan=64; %number of sites

%parameters for layer analysis
files(n).tip_loc_1=450;
files(n).tip_loc_2=425;
files(n).angle=45; 