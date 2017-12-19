ON=1; OFF=0;
%pathname = 'X:\Preycapture\';%miyazaki path to prey capture folder stored on Eyre 
%pathname = '/Users/jennifer/Dropbox/mucsimol_V1/Trial videos for Tristans tracking/Iphone5/'
%pathname = '/Users/jennifer/Desktop/green ear videos'
%pathname='H:\GreenEarTracking\iPhone\Good_video_tracks\Tracks\'
%pathname='S:\Tracks_11_22_17\'
pathname='C:\Users\Cris\Dropbox\Cell Types In Prey Capture\Tracks\Tracks_11_22_17\'
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
files(n).subj = '7587';%%bilateral hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Ntsr1_cre_CNO\Tracks_ALGS.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=1;
files(n).Moviefile='';

n=n+1;
files(n).subj = '7587';%%bilateral hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Ntsr1_cre_CNO\Tracks_ATPV.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=1;
files(n).Moviefile='';
% 
% n=n+1;
% files(n).subj = '7543';%video BMQM....
% files(n).sex = 'F';
% files(n).type=2;
% files(n).lighting = ON;
% files(n).Tfiles = 'Ntsr1_cre_CNO/Tracks_BMQM_virtual.mat';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=30;
% files(n).scale=22;%pixels/cm on video tracking
% files(n).group=1;
% files(n).Moviefile='';

n=n+1;
files(n).subj = 'NTSR1.3_LN';%bilateral hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Ntsr1_cre_CNO\Tracks_FNIB.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=1;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'NTSR1.3_LN';%bilateral hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Ntsr1_cre_CNO\Tracks_CCAH.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=1;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'NT-CH-1.2-LT';%bilateral hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Ntsr1_cre_CNO\Tracks_YAPF.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=1;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'NT-CH-1.2-LT';%bilateral hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Ntsr1_cre_CNO\Tracks_GBJV.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=1;
files(n).Moviefile='';


n=n+1;
files(n).subj = 'NT-Ch-1.6-LT';%bilateral hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Ntsr1_cre_CNO\Tracks_APBT_n';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=1;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'NT-Ch-1.6-LT';%bilateral hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Ntsr1_cre_CNO\Tracks_BXDA_n';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=1;
files(n).Moviefile='';

files(n).subj = 'NT-Ch-1.6-RT';%bilateral hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Ntsr1_cre_CNO\Tracks_IKUU_n';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=1;
files(n).Moviefile='';

files(n).subj = 'NT-Ch-1.6-RT';%bilateral hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Ntsr1_cre_CNO\Tracks_RYGH_n';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=1;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'NT-CH-1.6-LT';%bilateral hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Ntsr1_cre_CNO\Tracks_NXGU.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=1;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'NTSR1.1.7-LT';%bilateral hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Ntsr1_cre_CNO\Tracks_TXPP.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=1;
files(n).Moviefile='';
% 
n=n+1;
files(n).subj = 'NTSR1.1.8-LT';%bilateral hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Ntsr1_cre_CNO\Tracks_TYZB.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=1;
files(n).Moviefile='';
% 
n=n+1;
files(n).subj = 'NTSR1.1.8-LT';%bilateral hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Ntsr1_cre_CNO\Tracks_UWVJ.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=1;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'NTCH1.9-RT';%blateral hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Ntsr1_cre_CNO\Tracks_VZNJ.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=1;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'NTCH1.9-RT';%bilateral hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Ntsr1_cre_CNO\Tracks_YAWF.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=1;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GRP-1.7-LT';%bilateral hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GRP_cre_CNO\Tracks_OJIW.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=2;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GRP-1.7-LT';%bilateral hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GRP_cre_CNO\Tracks_DLWQ_g.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=2;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GRP_1.2_RT';
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GRP_cre_CNO\Tracks_XXMH.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=2;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GRP_1.2_RT';
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GRP_cre_CNO\Tracks_TLNS.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=2;
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
files(n).scale=22;%pixels/cm on video tracking
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
files(n).scale=22;%pixels/cm on video tracking
files(n).group=2;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GAi_1.4_RT';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GRP_cre_CNO\Tracks_HSTI.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=2;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GAi_1.4_RT';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GRP_cre_CNO\Tracks_SITB.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=2;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GAi_1.4_LN';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GRP_cre_CNO\Tracks_EMGR.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=2;
files(n).Moviefile='';


n=n+1;
files(n).subj = 'GAi_1.5_LT';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GRP_cre_CNO\Tracks_HRSX_g.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=2;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GAi_1.5_LT';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GRP_cre_CNO\Tracks_HTGT_g.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=2;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GAi_1.5_TT';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GRP_cre_CNO\Tracks_MAKZ_g.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=2;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GAi_1.5_TT';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GRP_cre_CNO\Tracks_QGRO_g.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=2;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GRP_1.4_LT';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GRP_cre_CNO\Tracks_TFEO_g.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=2;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GRP_1.4_LT';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GRP_cre_CNO\Tracks_TGQB_g.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=2;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GRP_1.4_TT';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GRP_cre_CNO\Tracks_TXGK_g.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=2;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GRP_1.4_TT';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GRP_cre_CNO\Tracks_XXVA_g.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=2;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GRP_1.5_LT';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GRP_cre_CNO\Tracks_TUMP.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=2;
files(n).Moviefile='';


n=n+1;
files(n).subj = 'PV_1.3_LT';
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'PV_cre_CNO\Tracks_CWWX.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
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
files(n).scale=22;%pixels/cm on video tracking
files(n).group=3;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'PV_1.2_TT';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'PV_cre_CNO\Tracks_QOZP.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=3;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'PV_1.2_TT';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'PV_cre_CNO\Tracks_QNZC.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=3;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'PV_1.4_LT';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'PV_cre_CNO\Tracks_JPAH.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=3;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'PV_1.4_LT';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'PV_cre_CNO\Tracks_IUSJ.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=3;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'PV_ChR_4B.17';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'PV_cre_CNO\Tracks_NOMG_p.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=3;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'PV_ChR_4B.17_e';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'PV_cre_CNO\Tracks_IFGS_p.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=3;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'PV_ChR_4B.17_e';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'PV_cre_CNO\Tracks_MLVA.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=3;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'PV_ChR_4B.17';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'PV_cre_CNO\Tracks_AHVM_p.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=3;
files(n).Moviefile='';
% 
n=n+1;
files(n).subj = 'PV_ChR_4B.17_2';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'PV_cre_CNO\Tracks_EOBT_p.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=3;
files(n).Moviefile='';
% 
n=n+1;
files(n).subj = 'PV_ChR_4B.17_2';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'PV_cre_CNO\Tracks_KFGM_p.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=3;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'PV_ChR_4B.17_A';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'PV_cre_CNO\Tracks_IWWU.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=3;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'PV_ChR_4B.17_A';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'PV_cre_CNO\Tracks_KWRS.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=3;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'PV_ChR_4B.17_B';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'PV_cre_CNO\Tracks_QJTK.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=3;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'PV_ChR_4B.17_B';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'PV_cre_CNO\Tracks_VRTJ.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=3;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'PV_ChR_4B.17_B';%hM4Di + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'PV_cre_CNO\Tracks_VWHD.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=3;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'Ntsr1_1.4_RT';%mcherry +CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\Tracks_BHDN.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'Ntsr1_1.4_RT'; %mcherry + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\Tracks_XKVV.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';


n=n+1;
files(n).subj = 'Grp 1.2-RT';%hM4_no expressoin+CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\Tracks_ENUW.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'Grp 1.2-RT';%hM4_no expressoin+CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\Tracks_EURG.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'PV_1.2_TT';%mcherry+CNO.
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\Tracks_CEXS.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'PV_1.2_TT';%mcherry+CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\Tracks_FDPE.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';
% 
n=n+1;
files(n).subj = 'GRP_1.4-LN';%mcherry + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\Tracks_FDPE.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GRP_1.4-LN';%mcherry + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\Tracks_AALS.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';


n=n+1;
files(n).subj = 'NTSR1.1.7-RT';%mcherry + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\Tracks_AMQU.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'NTSR1.1.7-RT';%mcherry + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\Tracks_CSJK.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';

n=n+1;
files(n).subj = '7583';%mcherry + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\BTQX_c.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';

n=n+1;
files(n).subj = '7583';%mcherry + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\MXHC_c.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';


n=n+1;
files(n).subj = '7581';%mcherry + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\Tracks_RHSO.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';


n=n+1;
files(n).subj = 'GRP 1.4-RN';%mcherry + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\Tracks_GFXD.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GRP 1.4-RN';%mcherry + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\Tracks_MXTH.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'PV-ChR2 4B.19-A';%mcherry + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\Tracks_QSBH.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'PV-ChR2 4B.19-A';%mcherry + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\Tracks_PCBC.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'PV-ChR2 4B.19-B';%mcherry + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'Control\POCZ_c.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=4;
files(n).Moviefile='';


n=n+1;
files(n).subj = 'GF1';%mcherry + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GF\Tracks_AGDW.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=5;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GF1';%mcherry + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GF\Tracks_FBYL.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=5;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GF1';%mcherry + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GF\Tracks_FEKB.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=5;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GF1';%mcherry + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GF\Tracks_HQHK.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=5;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GF1';%mcherry + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GF\Tracks_HQEE.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=5;
files(n).Moviefile='';

n=n+1;
files(n).subj = 'GF1';%mcherry + CNO
files(n).sex = 'F';
files(n).type=1;
files(n).lighting = ON;
files(n).Tfiles = 'GF\Tracks_MHJH.mat';
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=30;
files(n).scale=22;%pixels/cm on video tracking
files(n).group=5;
files(n).Moviefile='';

% n=n+1;
% files(n).subj = 'GF1';%mcherry + CNO
% files(n).sex = 'F';
% files(n).type=1;
% files(n).lighting = ON;
% files(n).Tfiles = 'GF\Tracks_MKKJ.mat';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=30;
% files(n).scale=22;%pixels/cm on video tracking
% files(n).group=5;
% files(n).Moviefile='';

% n=n+1;
% files(n).subj = 'GF1';%mcherry + CNO
% files(n).sex = 'F';
% files(n).type=1;
% files(n).lighting = ON;
% files(n).Tfiles = 'GF\Tracks_NLIG.mat';
% files(n).contrast = 100;
% files(n).notes = '';
% files(n).fps=30;
% files(n).scale=22;%pixels/cm on video tracking
% files(n).group=5;
% files(n).Moviefile='';
