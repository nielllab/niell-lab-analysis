% please provide only one session per datafile.
% this way we can organize things more smoothly.

global outputDir % this is where we will store the data, e.g. 'C:\\data' or ~/data
global inputDir
global info

outputDir = 'D:\Jen_analysis\NR5A_Pinping\Jen_NR5A_analysis_files\analysis_files\Martin_analysis';
inputDir  = 'D:/Jen_analysis/NR5A_Pinping/Jen_NR5A_analysis_files/analysis_files/Martin_analysis/Data/Hoy';

info      = [];

%%1==wt, 0==N2BKO, 2==N2AKO

% info(end+1).dataname = '2_2_15_analysis_2'; 
% info(end).datafile   = fullfile(inputDir, '2_2_15_analysis_2');
% info(end).stimDur    = 3.0535;
% info(end).stimType   = 'bars';
% info(end).genotype   = 0; % control;
% 
% info(end+1).dataname = '2_25_15_analysis_2'; 
% info(end).datafile   = fullfile(inputDir, '2_25_15_analysis_2');
% info(end).stimDur    = 3.0535;
% info(end).stimType   = 'bars';
% info(end).genotype   = 1;
% 
% info(end+1).dataname = '3_3_15_analysis_2'; 
% info(end).datafile   = fullfile(inputDir, '3_3_15_analysis_2');
% info(end).stimDur    = 3.0535;
% info(end).stimType   = 'bars';
% info(end).genotype   = 2;
% 
% info(end+1).dataname = '3_4_15_analysis_2'; 
% info(end).datafile   = fullfile(inputDir, '3_4_15_analysis_2');
% info(end).stimDur    = 3.0535;
% info(end).stimType   = 'bars';
% info(end).genotype   = 2;
% 
% info(end+1).dataname = '3_11_15_analysis_2'; 
% info(end).datafile   = fullfile(inputDir, '3_11_15_analysis_2');
% info(end).stimDur    = 3.0535;
% info(end).stimType   = 'bars';
% info(end).genotype   = 1;
% 
% info(end+1).dataname = '3_25_15_analysis_3_25_15'; 
% info(end).datafile   = fullfile(inputDir, '3_25_15_analysis_3_25_15');
% info(end).stimDur    = 3.0535;
% info(end).stimType   = 'bars';
% info(end).genotype   = 0;
% 
% info(end+1).dataname = '3_26_15_analysis_2'; 
% info(end).datafile   = fullfile(inputDir, '3_26_15_analysis_2');
% info(end).stimDur    = 3.0535;
% info(end).stimType   = 'bars';
% info(end).genotype   = 0;
% 
% info(end+1).dataname = '4_9_15_analysis_2'; 
% info(end).datafile   = fullfile(inputDir, '4_9_15_analysis_2');
% info(end).stimDur    = 3.0535;
% info(end).stimType   = 'bars';
% info(end).genotype   = 1;
% 
% info(end+1).dataname = '4_10_15_analysis_2'; 
% info(end).datafile   = fullfile(inputDir, '4_10_15_analysis_2');
% info(end).stimDur    = 3.0535;
% info(end).stimType   = 'bars';
% info(end).genotype   = 1;
% 
% info(end+1).dataname = '4_13_15_analysis_2'; 
% info(end).datafile   = fullfile(inputDir, '4_13_15_analysis_2');
% info(end).stimDur    = 3.0535;
% info(end).stimType   = 'bars';
% info(end).genotype   = 1;

info(end+1).dataname = '4_23_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '4_23_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

% info(end+1).dataname = '4_24_15_analysis_2A'; 
% info(end).datafile   = fullfile(inputDir, '4_24_15_analysis_2A');
% info(end).stimDur    = 3.0535;
% info(end).stimType   = 'bars';
% info(end).genotype   = 2;
% 
% info(end+1).dataname = '4_30_15_analysis_2'; 
% info(end).datafile   = fullfile(inputDir, '4_30_15_analysis_2');
% info(end).stimDur    = 3.0535;
% info(end).stimType   = 'bars';
% info(end).genotype   = 1;

info(end+1).dataname = '5_4_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '5_4_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 1;

info(end+1).dataname = '5_11_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '5_11_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 1;

info(end+1).dataname = '5_12_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '5_12_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;


info(end+1).dataname = '5_14_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '5_14_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 1;

info(end+1).dataname = '6_18_15_analysis_2A'; 
info(end).datafile   = fullfile(inputDir, '6_18_15_analysis_2A');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '6_22_15_analysis_2A'; 
info(end).datafile   = fullfile(inputDir, '6_22_15_analysis_2A');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '6_23_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '6_23_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '6_25_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '6_25_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 1;

info(end+1).dataname = '6_26_15_analysis_2A'; 
info(end).datafile   = fullfile(inputDir, '6_26_15_analysis_2A');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '6_28_15_analysis_2A'; 
info(end).datafile   = fullfile(inputDir, '6_28_15_analysis_2A');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;

info(end+1).dataname = '6_29_15_analysis_2A'; 
info(end).datafile   = fullfile(inputDir, '6_29_15_analysis_2A');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 1;

info(end+1).dataname = '6_30_15_analysis_2B'; 
info(end).datafile   = fullfile(inputDir, '6_30_15_analysis_2B');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;

info(end+1).dataname = '7_1_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '7_1_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '7_29_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '7_29_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;

info(end+1).dataname = '8_7_15_analysis_2A'; 
info(end).datafile   = fullfile(inputDir, '8_7_15_analysis_2A');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;

info(end+1).dataname = '8_10_15_analysis_2A'; 
info(end).datafile   = fullfile(inputDir, '8_10_15_analysis_2A');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 1;

info(end+1).dataname = '8_11_15_analysis_2A'; 
info(end).datafile   = fullfile(inputDir, '8_11_15_analysis_2A');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '8_11_15_rec2_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '8_11_15_rec2_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '8_13_15_analysis_2A'; 
info(end).datafile   = fullfile(inputDir, '8_13_15_analysis_2A');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '8_14_15_analysis'; 
info(end).datafile   = fullfile(inputDir, '8_14_15_analysis');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '8_17_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '8_17_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 1;

info(end+1).dataname = '8_18_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '8_18_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;

info(end+1).dataname = '8_18_15_rec2_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '8_18_15_rec2_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;

info(end+1).dataname = '8_19_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '8_19_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;


info(end+1).dataname = '8_19_15_rec2_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '8_19_15_rec2_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 0;

info(end+1).dataname = '9_12_15_analysis'; 
info(end).datafile   = fullfile(inputDir, '9_12_15_analysis');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;

info(end+1).dataname = '9_16_15_analysis_2'; 
info(end).datafile   = fullfile(inputDir, '9_16_15_analysis_2');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;

info(end+1).dataname = 'analysis_9_14_15'; 
info(end).datafile   = fullfile(inputDir, 'analysis_9_14_15');
info(end).stimDur    = 1.5;
info(end).stimType   = 'drift';
info(end).genotype   = 2;
