% please provide only one session per datafile.
% this way we can organize things more smoothly.

global outputDir % this is where we will store the data, e.g. 'C:\\data' or ~/data
global inputDir
global info

outputDir = 'F:\Jennifer_Development\Jen_analysis_development\Data_2_17_17\data';
inputDir  = 'F:/Jennifer_Development/Jen_analysis_development/Data_2_17_17/';

info      = [];

%%1==adult, 2==EO1

info(end+1).dataname = 'analysis_5_25_13_B_rec1_adult'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_5_25_13_B_rec1_adult');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1; % Adult;
% 
info(end+1).dataname = 'analysis_7_25_13_mouseB_rec1'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_7_25_13_mouseB_rec1');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;
% 
info(end+1).dataname = 'analysis_8_7_13_EO1_rec1_more_strigent'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_8_7_13_EO1_rec1_more_strigent');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;
% 
info(end+1).dataname = 'analysis_8_7_13_rec2'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_8_7_13_rec2');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;
% 
info(end+1).dataname = 'analysis_8_8_13_rec1'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_8_8_13_rec1');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;
% 
info(end+1).dataname = 'analysis_9_9_13_EO1_rec1'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_9_9_13_EO1_rec1');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;
% 
info(end+1).dataname = 'analysis_9_9_13_EO1_rec2'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_9_9_13_EO1_rec2');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;
% 
info(end+1).dataname = 'analysis_9_12_13_adult_rec1'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_9_12_13_adult_rec1');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;
% 
info(end+1).dataname = 'analysis_9_30_13_EO1_rec1'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_9_30_13_EO1_rec1');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;
% 
info(end+1).dataname = 'analysis_9_30_13_EO1_rec2'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_9_30_13_EO1_rec2');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;

info(end+1).dataname = 'analysis_11_11_13_adult_rec1'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_11_11_13_adult_rec1');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;

info(end+1).dataname = 'analysis_11_11_13_adult_rec2'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_11_11_13_adult_rec2');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;

info(end+1).dataname = 'analysis_11_13_13_rec2'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_11_13_13_rec2');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;

info(end+1).dataname = 'analysis_11_14_13_adult_rec1'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_11_14_13_adult_rec1');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;

info(end+1).dataname = 'analysis_11_14_13_adult_rec2'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_11_14_13_adult_rec2');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;

info(end+1).dataname = 'analysis_11_15_13_adult_rec1'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_11_15_13_adult_rec1');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;


info(end+1).dataname = 'analysis_11_15_13_adult_rec2'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_11_15_13_adult_rec2');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;

info(end+1).dataname = 'analysis_11_20_13_EO0_rec2'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_11_20_13_EO0_rec2');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;

info(end+1).dataname = 'analysis_11_20_13_rec1'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_11_20_13_rec1');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;

info(end+1).dataname = 'analysis_11_21_13_mouseA_EO1'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_11_21_13_mouseA_EO1');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;

info(end+1).dataname = 'analysis_11_23_13_EO2'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_11_23_13_EO2');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;

info(end+1).dataname = 'analysis_11_26_13_EO1_mouseA_rec1'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_11_26_13_EO1_mouseA_rec1');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;

info(end+1).dataname = 'analysis_11_26_13_EO1_rec2'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_11_26_13_EO1_rec2');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;

info(end+1).dataname = 'analysis_11_26_13_mouseB_EO1'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_11_26_13_mouseB_EO1');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;

info(end+1).dataname = 'analysis_adult_11_13_13_rec1'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_adult_11_13_13_rec1');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;

info(end+1).dataname = 'analysis_cluster_data_07_25_13_mouseB_adult_rec2'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_cluster_data_07_25_13_mouseB_adult_rec2');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;

info(end+1).dataname = 'analysis_rec1_A_5_22_13_strict_selection'; 
info(end).datafile   = fullfile(inputDir, 'Bar/analysis_rec1_A_5_22_13_strict_selection');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 2;

