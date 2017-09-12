ON=1; OFF=0;
pathname = 'W:\';

%miyazaki path to prey capture folder stored on Eyre
%pathname='E:\Prey_capture_data\'; %eyre path

% if ismac
%     pathname = '/Users/crisniell/Dropbox/Prey capture/Prey_tracks/'
% end
% 
% pathname = '/Users/crisniell/Dropbox/Prey capture/Prey_tracks/'
% 

n=0;


%%%group: 1=wt blk6; 2= Dark 1st exp; 3=eyeSuture; 4=Bbgd; 5=earPlug;
% 6=EyeReopen; 7=EarUnplugged; 8=blind mice; 10=ctl Blind
%%%female=1 and male=2
box_Bounds(:,1)=[231.546931	1060.196751	1676.976534	1063.893502	238.940433 6.622744 1669.583032	14.016245];
box_Bounds(:,1)=[283.301444	1067.590253	312.875451	17.712996	1721.337545	1074.983755	1725.034296	10.319495];

box_Bounds(:,2)=[198.276173	1089.770758	1850.723827	1078.680505	194.579422 0 1861.814079 2.925993];
box_Bounds(:,3)=[283.30144 1086.074007 1832.240072 1086.074007 279.604693 0 1850.723827	0];
 figure;plot(box_Bounds([1 3 7 5 1],1),box_Bounds([2 4 8 6 2],1));% plots
%  boundry #1


%%

% approach to virtual stim

n=n+1;
files(n).subj = 'Prey24_NT';
files(n).lighting = ON;
files(n).trackpts = 'Size_9_1_16\Prey24_NT_slide11\DLTdv5_data_xypts.csv';%specific path and file name
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;
files(n).FrameS=535;%frame where cricket is first available
files(n).FrameEnd=740;%frame where cricket is caught
files(n).CapTime=11;
files(n).FrameFlash='';
files(n).Moviefile='Size_9_1_16\Prey24_NT_slide11.avi';
files(n).body=0;
files(n).bounds=box_Bounds(:,1);
%%


