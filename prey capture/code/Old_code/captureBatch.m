ON=1; OFF=0;
pathname = 'F:Prey_tracks\';
if ismac
    pathname = '/Users/crisniell/Dropbox/Prey capture/Prey_tracks'
end


n=0;

%%%% prey 1
%%%group: 1=wt blk6 mice; 2= EyeSuture; 3=BlkBgd; 4=EarPlug; 5=trimmed whiskers;
%%%6=Plexiglass; 7=EyeReopen; 8=EarUnplugged; 9=blind mice; 10=ctl Blind

n=n+1;
files(n).subj = 'prey1';
files(n).lighting = OFF;
files(n).trackpts = 'Prey1\7_28_15\Prey1_7_28_15_trial2_NIR.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey1';
files(n).lighting = ON;
files(n).trackpts = 'Prey1\7_28_15\prey1_7_28_15_trial3_FS.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey1';
files(n).lighting = OFF;
files(n).trackpts = 'Prey1\7_28_15\Prey1_7_28_15_trial4_NIR.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey1';
files(n).lighting = ON;
files(n).trackpts = 'Prey1\7_28_15\Prey1_7_28_15_trial5_FS.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey1';
files(n).lighting = ON;
files(n).trackpts = 'Prey1\9_2_15\Prey1_9_2_15_trial1_FS.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey1';
files(n).lighting = OFF;
files(n).trackpts = 'Prey1\9_2_15\Prey1_9_2_15_trial2_NIR.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;


%%% prey 2

n=n+1;
files(n).subj = 'prey2';
files(n).lighting = ON;
files(n).trackpts = 'Prey2\9_16_15\9_16_15_trial2_FS.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey2';
files(n).lighting = OFF;
files(n).trackpts = 'Prey2\9_16_15\9_16_15_trial3_NIR.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey2';
files(n).lighting = ON;
files(n).trackpts = 'Prey2\9_16_15\9_16_15_trial6_FS.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey2';
files(n).lighting = OFF;
files(n).trackpts = 'Prey2\9_16_15\9_16_15_trial7_NIR.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey2';
files(n).lighting = OFF;
files(n).trackpts = 'Prey2\9_16_15\9_16_15_trial8_NIR.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey2';
files(n).lighting = ON;
files(n).trackpts = 'Prey2\9_16_15\9_16_15_trial9_FS.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=1;

%%%Prey4_LT

n=n+1;
files(n).subj = 'prey4_LT';
files(n).lighting = OFF;
files(n).trackpts = 'Prey4_LT\11_6_15\Prey4_LT_11_6_15_trial4_NIR.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4_LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_LT\11_6_15\Prey4_LT_11_6_15_trial5_FS.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4_LT-F';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_LT\11_6_15\Prey4_LT_11_6_15_trial7_FS.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4_LT-N';
files(n).lighting = OFF;
files(n).trackpts = 'Prey4_LT\11_6_15\Prey4_LT_11_6_15_trial8_NIR.csv';
files(n).soundfile = '';
files(n).headbar = 1;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4_RT-N';
files(n).lighting = OFF;
files(n).trackpts = 'Prey4_RT\11_6_15\Prey4_RT_11_6_15_trial2_NIR.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4_RT-F';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_RT\11_6_15\Prey4_RT_11_6_15_trial3_FS.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4_RT-N';
files(n).lighting = OFF;
files(n).trackpts = 'Prey4_RT\11_6_15\Prey4_RT_11_6_15_trial4_NIR.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4_RT-F';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_RT\11_6_15\Prey4_RT_11_6_15_trial5_FS.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=2;

%Prey4_LN

n=n+1;
files(n).subj = 'prey4_LN-N';
files(n).lighting = OFF;
files(n).trackpts = 'Prey4_LN\11_7_15\11_7_15_Prey4-lN_Trial2_NIR.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4_LN-F';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_LN\11_7_15\11_7_15_Prey4_LN_Trial5_FS.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4_LN-N';
files(n).lighting = OFF;
files(n).trackpts = 'Prey4_LN\11_7_15\11_7_15_Prey4_LN_Trial4_NIR.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=1;
files(n).sex=2;

%prey6 blind mice
n=n+1;
files(n).subj = 'prey6_RT_FS';
files(n).lighting = ON;
files(n).trackpts = 'Blind_prey6_RT\track_1_15_16_prey6_RT_trial1_FS.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=9;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey6_RT_NIR';
files(n).lighting = OFF;
files(n).trackpts = 'Blind_prey6_RT\track_1_15_16_prey6_RT_trial2_NIR.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=9;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey6_RT_FS';
files(n).lighting = ON;
files(n).trackpts = 'Blind_prey6_RT\track_1_15_16_prey6_RT_trial3_FS.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=9;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey6_RT_NIR';
files(n).lighting = OFF;
files(n).trackpts = 'Blind_prey6_RT\track_1_15_16_prey6_RT_trial4_NIR.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=9;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey6_RT_FS';
files(n).lighting = ON;
files(n).trackpts = 'Blind_prey6_RT\track_1_15_16_prey6_RT_trial5_FS.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=9;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey6_LT_FS';
files(n).lighting = ON;
files(n).trackpts = 'Blind_prey6_LT\track_1_15_16_prey6_LT_trial1_FS.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=9;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey6_LT_NIR';
files(n).lighting = ON;
files(n).trackpts = 'Blind_prey6_LT\track_1_15_16_prey6_LT_trial2_NIR.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=9;
files(n).sex=2;

%% BBgd experiments
n=n+1;
files(n).subj = 'prey4-LN';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_LN\BBgd\Prey4_LN_BBgd_FS_trial5_2_5_16.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=3;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_LT\BBgd\Prey4_LT_BBgd_FS_trial4_2_5_16.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=3;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_LT\BBgd\Prey4_LT_BBgd_FS_trial6_2_5_16.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=3;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-RT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_RT\BBgd\Prey4_RT_BBgd_FS_trial3_2_5_16.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=3;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-RT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_RT\BBgd\Prey4_RT_BBgd_FS_trial4_2_5_16.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=3;
files(n).sex=2;


%%%group: 1=wt blk6 mice; 2= EyeSuture; 3=BlkBgd; 4=EarPlug; 5=trimmed whiskers;
%%%6=Plexiglass; 7=EyeReopen; 8=EarUnplugged; 9=blind mice; 10=ctl Blind

%% Eye Suture experiments

n=n+1;
files(n).subj = 'prey4-LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_LT\EyeSuture\Prey4_LT_ES_FS_trial2_2_8_16.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=2;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_LT\EyeSuture\Prey4_LT_ES_FS_trial4_2_8_16.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=2;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_LT\EyeSuture\Prey4_LT_FS_ES_trial7.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=2;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_LT\EyeSuture\Prey4_LT_ES_FS_trial5_2_8_16.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=2;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_LT\EyeSuture\Prey4_LT_ES_FS_trial6_2_8_16.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=2;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-LT';
files(n).lighting = OFF;
files(n).trackpts = 'Prey4_LT\EyeSuture\Prey4_LT_ES_NIR_trial3_2_8_16.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=2;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-TT';
files(n).lighting = OFF;
files(n).trackpts = 'Prey4_TT\EyeSuture\Prey4_TT_NIR_ES_trial1_noBody.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=2;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-TT';
files(n).lighting = OFF;
files(n).trackpts = 'Prey4_TT\EyeSuture\Prey4_TT_NIR_ES_trial4_noBody.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=2;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-TT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_TT\EyeSuture\Prey4_TT_FS_ES_trial1.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=2;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey5-LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey5_LT\EyeSuture\Prey5_LT_ES_FS_trial1.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=2;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey5-LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey5_LT\EyeSuture\Prey5_LT_ES_FS_trial2.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=2;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey5-LT';
files(n).lighting = OFF;
files(n).trackpts = 'Prey5_LT\EyeSuture\Prey5_LT_ES_NIR_trial1.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=2;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey5-LT';
files(n).lighting = OFF;
files(n).trackpts = 'Prey5_LT\EyeSuture\Prey5_LT_ES_NIR_trial3.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=2;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey5-LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey5_LT\EyeSuture\Prey5_LT_ES_NIR_trial4.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=2;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey5-LT';
files(n).lighting = OFF;
files(n).trackpts = 'Prey5_LT\EyeSuture\Prey5_LT_ES_NIR_trial5.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=2;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey5-LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey5_LT\EyeSuture\Prey5_LT_ES_FS_trial6.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 50;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=2;
files(n).sex=1;

% Ear Plugged

n=n+1;
files(n).subj = 'prey4-LN';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_LN\EarPlug\Prey4_LN_FS_EP_trial1.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-LN';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_LN\EarPlug\Prey4_LN_FS_EP_trial3_data_xypts.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-LN';
files(n).lighting = OFF;
files(n).trackpts = 'Prey4_LN\EarPlug\Prey4_LN_NIR_EP_trial2.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-LN';
files(n).lighting = OFF;
files(n).trackpts = 'Prey4_LN\EarPlug\Prey4_LN_NIR_EP_trial4.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-LT';
files(n).lighting = OFF;
files(n).trackpts = 'Prey4_LT\EarPlug\Prey4_LT_NIR_EP_trial7_noBody.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_LT\EarPlug\Prey4_LT_FS_EP_trial9.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_LT\EarPlug\Prey4_LT_FS_EP_trial10.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-RT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_RT\EarPlug\Prey4_RT_FS_EP_trial1_noBody.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-RT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_RT\EarPlug\Prey4_RT_FS_EP_trial2_noBody.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-RT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_RT\EarPlug\Prey4_RT_FS_EP_trial3_noBody.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-RT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_RT\EarPlug\Prey4_RT_FS_EP_trial4_noBody.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-RT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_RT\EarPlug\Prey4_RT_FS_EP_trial5_noBody.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-RT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_RT\EarPlug\Prey4_RT_FS_EP_trial6_noBody.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-RT';
files(n).lighting = OFF;
files(n).trackpts = 'Prey4_RT\EarPlug\Prey4_RT_NIR_EP_trial6_noBody.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey5-LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey5_LT\EarPlug\Prey5_LT_FS_EP_trial1.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey5-LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey5_LT\EarPlug\Prey5_LT_FS_EP_trial3.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey5-LT';
files(n).lighting = OFF;
files(n).trackpts = 'Prey5_LT\EarPlug\Prey5_LT_NIR_EP_trial2.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey5-LT';
files(n).lighting = OFF;
files(n).trackpts = 'Prey5_LT\EarPlug\Prey5_LT_NIR_EP_trial4.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey5-TT';
files(n).lighting = ON;
files(n).trackpts = 'Prey5_TT\Earplug\Prey5_TT_FS_EP_trial2_noBody.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey5-TT';
files(n).lighting = ON;
files(n).trackpts = 'Prey5_TT\Earplug\Prey5_TT_FS_EP_trial6_noBody.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey5-TT';
files(n).lighting = OFF;
files(n).trackpts = 'Prey5_TT\Earplug\Prey5_TT_NIR_trial4_noBody.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=4;
files(n).sex=1;


%% eye suture removal
n=n+1;
files(n).subj = 'prey4-LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_LT\ESremoved\Prey4_LT_ESR_FS_trial1.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=7;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_LT\ESremoved\Prey4_LT_ESR_FS_trial2.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=7;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey4_LT\ESremoved\Prey4_LT_ESR_FS_trial3.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=7;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-LT';
files(n).lighting = OFF;
files(n).trackpts = 'Prey4_LT\ESremoved\Prey4_LT_ESR_NIR_trial2.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=7;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey4-LT';
files(n).lighting = OFF;
files(n).trackpts = 'Prey4_LT\ESremoved\Prey4_LT_ESR_NIR_trial4.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=7;
files(n).sex=2;

n=n+1;
files(n).subj = 'prey5-LT';
files(n).lighting = OFF;
files(n).trackpts = 'Prey5_LT\EyeSutureRemoved\Prey5_LT_NIR_ESR_trial1_noBody.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=7;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey5-LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey5_LT\EyeSutureRemoved\Prey5_LT_FS_ESR_trial1_noBody.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=7;
files(n).sex=1;

n=n+1;
files(n).subj = 'prey5-LT';
files(n).lighting = ON;
files(n).trackpts = 'Prey5_LT\EyeSutureRemoved\Prey5_LT_FS_ESR_trial3_noBody.csv';
files(n).soundfile = '';
files(n).headbar = 0;
files(n).contrast = 100;
files(n).notes = '';
files(n).fps=60;
files(n).scale=35;%pixels/cm on video tracking
files(n).group=7;
files(n).sex=1;

%%%group: 1=wt blk6 mice; 2= EyeSuture; 3=BlkBgd; 4=EarPlug; 5=trimmed whiskers;
%%%6=Plexiglass; 7=EyeReopen; 8=EarUnplugged; 9=blind mice; 10=ctl Blind

%%no whiskers
% n=n+1;
% files(n).subj = 'prey1-F-NW';
% files(n).lighting = ON;
% files(n).trackpts = 'Prey1\11_11_15_nowhisk\prey1_11_11_15_nowhisk_trial1_FS.csv';
% files(n).soundfile = '';
% files(n).headbar = 0;
% files(n).contrast = 100;
% files(n).notes = 'nowhisk';
% files(n).fps=60;
% 





