%prey capture recipees, work flow
%% generate track data using DLTdv5 from Ty Hedricks group. code can be downloaded and is maintained on Ty's website.

DLTdv5.m

%%create a "batch file" where all input variables for each file are organized and entered. For example,scale of arena, subject ID, file path and number of frames analyzed. X = good name for your overall experiment experiment 

captureBatchX.m

%% run compile script to call subscripts that generate relevant measures such as time to capture, head angle relative to prey etc. And then performs plotting functions used in the original paper

CompileCaptureLoopsX.m

%% functions and outside scripts called within compile code
analyzeCapture.m
getSmoothAngle.m
loadMouseMovie.m
trackMovie.m
myHist2.m


