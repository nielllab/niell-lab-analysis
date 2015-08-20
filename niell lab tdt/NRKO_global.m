% please provide only one session per datafile.
% this way we can organize things more smoothly.

global outputDir % this is where we will store the data, e.g. 'C:\\data' or ~/data
global inputDir
global info

outputDir = 'D:/Jen_analysis/analysis_martin/Data/Hoy/Analysis';
inputDir  = 'D:/Jen_analysis/analysis_martin/Data/Hoy';

info      = [];
info(end+1).dataname = '2_19_15'; 
info(end).datafile   = fullfile(inputDir, '2_19_15');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 0; % control;

info(end+1).dataname = '2_23_15_JLH'; 
info(end).datafile   = fullfile(inputDir, '2_23_15_JLH');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;

info(end+1).dataname = '2_24_15'; 
info(end).datafile   = fullfile(inputDir, '2_24_15');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;

info(end+1).dataname = '2_25_15'; 
info(end).datafile   = fullfile(inputDir, '2_25_15');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 0;

info(end+1).dataname = '2_2_15'; 
info(end).datafile   = fullfile(inputDir, '2_2_15');
info(end).stimDur    = 3.0535;
info(end).stimType   = 'bars';
info(end).genotype   = 1;
