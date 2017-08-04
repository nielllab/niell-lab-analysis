ON=1; OFF=0;
%pathname = 'X:\Preycapture\';%miyazaki path to prey capture folder stored on Eyre 
%pathname = '/Users/jennifer/Dropbox/mucsimol_V1/Trial videos for Tristans tracking/Iphone5/'
%pathname = '/Users/jennifer/Desktop/green ear videos'
%pathname='F:\GreenEarTracking\iPhone\Good_video_tracks\Tracks\';%on Jen-W10
pathname='Y:\GreenEarTracking\iPhone\Good_video_tracks\Tracks\';%on Niell-V2

n=0;


%%%female=1 and male=2;
%%%contrast:100
%%%size:100
%%%age: adult=1, juvenile= 2-->5
%%%grouplabels = {'Ntsr1_hM4+CNO','GRP_hM4+CNO','PV_hM4+CNO','cre_mCherry+CNO','cre_hM4_saline'};

%%

%pathname='F:\GreenEarTracking\iPhone'

%type 1= live type 2 = virtual

n=n+1;
files(n).subj = '7587';%video ALGS010....
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Ntsr1_cre_CNO\Tracks_ALGS.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=21;%pixels/cm on video tracking
files(n).group=1;
files(n).Moviefile='';
% 
n=n+1;
files(n).subj = '7587';%video ALGS010....
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Ntsr1_cre_CNO\Tracks_NXGU.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=21;%pixels/cm on video tracking
files(n).group=1;
files(n).Moviefile='';
% n=n+1;
% files(n).subj = '7543';%video BMQM....
% files(n).sex = 'F';
% files(n).type=2;
% files(n).lighting = ON;
% files(n).Tfiles = 'Ntsr1_cre_CNO/Tracks_BMQM_virtual.mat';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=30;
% files(n).scale=21;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='';

n=n+1;
files(n).subj = 'NTSR1.3_LN';
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Ntsr1_cre_CNO\Tracks_FNIB.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=21;%pixels/cm on video tracking
files(n).group=1;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'NCH-1.6_LT';
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Ntsr1_cre_CNO\Tracks_GBJV.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=21;%pixels/cm on video tracking
files(n).group=1;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GRP_1.2_LT';
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GRP_cre_CNO\Tracks_QPDL.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=21;%pixels/cm on video tracking
files(n).group=2;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GRP_1.2_LT';
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GRP_cre_CNO\Tracks_MISG.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=21;%pixels/cm on video tracking
files(n).group=2;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GAi-RT';
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GRP_cre_CNO\Tracks_HSTI.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=21;%pixels/cm on video tracking
files(n).group=2;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GAi-RT';
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GRP_cre_CNO\Tracks_TLNS.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=21;%pixels/cm on video tracking
files(n).group=2;
files(n).Moviefile='';


n=n+1;
files(n).subj = 'PV_1.2_TT';%CEXS.
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'PV_cre_CNO\Tracks_CEXS.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=21;%pixels/cm on video tracking
files(n).group=3;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'PV_1.3_RT';%
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'PV_cre_CNO\Tracks_CWWX.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=21;%pixels/cm on video tracking
files(n).group=3;
files(n).Moviefile='';


n=n+1;
files(n).subj = 'PV_1.3_LT';
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'PV_cre_CNO\Tracks_TLMM.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=21;%pixels/cm on video tracking
files(n).group=3;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'PV_1.3_LT';
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'PV_cre_CNO\Tracks_QOZP.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=21;%pixels/cm on video tracking
files(n).group=3;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'Ntsr1_1.4_RT';
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\Tracks_BHDN.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=21;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'Ntsr1_1.4_RT';
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\Tracks_XKVV.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=21;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GRP-RT_noExp';
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\Tracks_ENUW.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=21;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GRP-RT_noExp';
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\Tracks_EURG.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=21;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';







