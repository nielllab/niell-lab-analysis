% please provide only one session per datafile.
% this way we can organize things more smoothly.

global outputDir % this is where we will store the data, e.g. 'C:\\data' or ~/data
global inputDir
global info

outputDir = 'D:\Jen_analysis\analysis_martin\data_1_25_16\analysis_files';
inputDir  = 'D:/Jen_analysis/analysis_martin/data_1_25_16';

info      = [];

%% 0==N2BKO,1==wt, 2==N2AKO
info(end+1).dataname = '1_18_16_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/1_18_16_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;

info(end+1).dataname = '2_2_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/2_2_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 0;
% 
info(end+1).dataname = '2_25_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/2_25_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;
% 
info(end+1).dataname = '3_3_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/3_3_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;
% 
info(end+1).dataname = '3_4_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/3_4_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;
% 
info(end+1).dataname = '3_11_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/3_11_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;
% 
info(end+1).dataname = '3_25_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/3_25_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 0;
% 
info(end+1).dataname = '3_26_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/3_26_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 0;
% 
info(end+1).dataname = '4_9_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/4_9_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;
% 
info(end+1).dataname = '4_10_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/4_10_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;
% 
info(end+1).dataname = '4_13_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/4_13_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;

info(end+1).dataname = '4_23_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/4_23_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 0;

info(end+1).dataname = '4_24_15_analysis_2A_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/4_24_15_analysis_2A_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;

info(end+1).dataname = '4_30_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/4_30_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;

info(end+1).dataname = '5_4_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/5_4_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;

info(end+1).dataname = '5_11_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/5_11_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;

info(end+1).dataname = '5_12_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/5_12_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;


% info(end+1).dataname = '5_14_15_analysis_2_bar'; 
% info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/5_14_15_analysis_2_bar');
% info(end).stimDur    = 3.0535;
% info(end).stimType   = 'bars';
% info(end).genotype   = 1;

info(end+1).dataname = '6_18_15_analysis_2A_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/6_18_15_analysis_2A_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 0;

info(end+1).dataname = '6_22_15_analysis_2A_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/6_22_15_analysis_2A_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 0;

info(end+1).dataname = '6_23_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/6_23_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 0;

info(end+1).dataname = '6_25_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/6_25_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;

info(end+1).dataname = '6_26_15_analysis_2A_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/6_26_15_analysis_2A_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 0;

info(end+1).dataname = '6_28_15_analysis_2A_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/6_28_15_analysis_2A_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;

info(end+1).dataname = '6_29_15_analysis_2A_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/6_29_15_analysis_2A_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;

info(end+1).dataname = '6_30_15_analysis_2B_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/6_30_15_analysis_2B_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;

info(end+1).dataname = '7_1_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/7_1_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 0;

info(end+1).dataname = '7_29_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/7_29_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;

info(end+1).dataname = '8_11_15_analysis_2A_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/8_11_15_analysis_2A_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 0;


info(end+1).dataname = '8_13_15_analysis_2A_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/8_13_15_analysis_2A_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 0;

info(end+1).dataname = '8_14_15_analysis_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/8_14_15_analysis_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 0;

info(end+1).dataname = '8_17_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/8_17_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;

info(end+1).dataname = '8_18_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/8_18_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;

info(end+1).dataname = '8_18_15_rec2_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/8_18_15_rec2_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;

info(end+1).dataname = '8_19_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/8_19_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 0;


info(end+1).dataname = '8_19_15_rec2_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/8_19_15_rec2_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 0;

info(end+1).dataname = '9_12_15_analysis_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/9_12_15_analysis_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;

info(end+1).dataname = '9_16_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/9_16_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;

info(end+1).dataname = '11_2_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/11_2_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;

info(end+1).dataname = '11_3_15_analysis_2_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/11_3_15_analysis_2_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;

info(end+1).dataname = 'analysis_12_16_15_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/analysis_12_16_15_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 0;

info(end+1).dataname = 'analysis_1_20_16_bar'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_bar_data/analysis_1_20_16_bar');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;

%%begin drift
info(end+1).dataname = '1_13_15_analysis_pos_ctl'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/1_13_15_analysis_pos_ctl');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 1

info(end+1).dataname = '1_18_16_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/1_18_16_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 1

info(end+1).dataname = '2_2_15_analysis_2.mat'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/2_2_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0

info(end+1).dataname = '2_25_15_analysis_2.mat'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/2_25_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 1

info(end+1).dataname = '3_3_15_analysis_2.mat'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/3_3_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2

info(end+1).dataname = '3_4_15_analysis_2.mat'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/3_4_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2

info(end+1).dataname = '3_11_15_analysis_2.mat'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/3_11_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 1

info(end+1).dataname = '3_25_15_analysis_3_25_15.mat'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/3_25_15_analysis_3_25_15');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0

info(end+1).dataname = '3_26_15_analysis_2.mat'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/3_26_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0

info(end+1).dataname = '4_9_15_analysis_2.mat'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/4_9_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 1

info(end+1).dataname = '4_10_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/4_10_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 1

info(end+1).dataname = '4_13_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/4_13_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 1;

info(end+1).dataname = '4_23_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/4_23_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '4_24_15_analysis_2A'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/4_24_15_analysis_2A');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;

info(end+1).dataname = '4_30_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/4_30_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 1;

info(end+1).dataname = '5_4_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/5_4_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 1;

info(end+1).dataname = '5_11_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/5_11_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 1;

info(end+1).dataname = '5_12_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/5_12_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;


info(end+1).dataname = '6_18_15_analysis_2A'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/6_18_15_analysis_2A');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '6_22_15_analysis_2A'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/6_22_15_analysis_2A');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '6_23_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/6_23_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '6_25_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/6_25_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 1;

info(end+1).dataname = '6_26_15_analysis_2A'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/6_26_15_analysis_2A');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '6_28_15_analysis_2A'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/6_28_15_analysis_2A');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;

info(end+1).dataname = '6_29_15_analysis_2A'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/6_29_15_analysis_2A');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 1;

info(end+1).dataname = '6_30_15_analysis_2B'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/6_30_15_analysis_2B');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;

info(end+1).dataname = '7_1_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/7_1_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '7_29_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/7_29_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;

info(end+1).dataname = '8_11_15_analysis_2A'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/8_11_15_analysis_2A');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '8_11_15_rec2_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/8_11_15_rec2_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '8_13_15_rec2_analysis_2A'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/8_13_15_rec2_analysis_2A');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '8_13_15_analysis_2A'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/8_13_15_analysis_2A');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '8_14_15_analysis'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/8_14_15_analysis');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '8_17_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/8_17_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 1;

info(end+1).dataname = '8_18_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/8_18_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;

info(end+1).dataname = '8_18_15_rec2_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/8_18_15_rec2_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;

info(end+1).dataname = '8_19_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/8_19_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;


info(end+1).dataname = '8_19_15_rec2_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/8_19_15_rec2_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '9_12_15_analysis'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/9_12_15_analysis');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;

info(end+1).dataname = '9_16_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/9_16_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;

info(end+1).dataname = '11_2_5_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/11_2_5_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;

info(end+1).dataname = '11_3_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/11_3_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;

info(end+1).dataname = 'analysis_12_16_15'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/analysis_12_16_15');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = 'analysis_1_20_16'; 
info(end).datafile   = fullfile(inputDir, '/Hoy_drift_data/analysis_1_20_16');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;