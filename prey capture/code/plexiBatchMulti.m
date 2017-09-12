ON=1; OFF=0;
pathname = 'W:\';%miyazaki path to prey capture folder stored on Eyre 
n=0;


%%%group: 1=; 2= ; 3=; 4=; 5=;
%%%female=1 and male=2;
%%%contrast:100,50, 25, 12.5, and 6.25, percent transmittance
%%%size:100, 50, 25, 12.5, 6.25, scaled down size along all dimensions of
%%%ellipse (naturalistic stimuli)

n=n+1;
files(n).subj = '';
files(n).lighting = ON;
files(n).trackpts = '\DLTdv5_data_xypts.csv';%specific path and file name
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100; %enter grey filter ID
files(n).size = 100; % enter scaled virtual stimuli ID
files(n).notes = '';
files(n).fps=60; 
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1; %assign difference grey values to group numbers 
files(n).sex=1;
files(n).FrameS=407; %frame where cricket is first available
files(n).FrameEnd=1079; %frame where cricket is caught
files(n).CapTime=11;
files(n).FrameFlash='';
files(n).Moviefile='.avi';
files(n).body=0;
% 