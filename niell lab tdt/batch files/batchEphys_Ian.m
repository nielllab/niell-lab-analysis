pathname = '\\niell-v2-w7\Ian_analysis\';
n=0;

%% example entry%%
n=n+1;
files(n).expt = '071416';
files(n).dir = '07_14_16';
files(n).tank = '';
files(n).clusterfile = 'cluster_data_07_14_16_cluster_071416';
files(n).analysisfile = 'analysis_071416';
files(n).predark = 'dark_pre1';
files(n).postdark = 'dark_post1';
files(n).prewn = 'wn_pre1';
files(n).postwn = 'wn_post1';  
files(n).predrift = 'drift_pre1';
files(n).postdrift = 'drift_post1';
files(n).prebars = '';
files(n).postbars = '';
files(n).pretreatment = '';
files(n).injection = 'injection';
files(n).treatment = 'MGluR2'; 
files(n).prepinpFile = ''
files(n).postpinpFile = ''
files(n).layers = [];
files(n).notes = 'good data';
files(n).misc = '15mg/kg ly379268'; %dose
files(n).tip1 = 700; files(n).tip2 = 750; files(n).angle = 45;  %%% example
% files(n).predark_camera = 'dark_pre_051216_eye';  
% files(n).postdark_camera = 'dark_post_051216_eye';
files(n).badsites = [];




for f = 1:length(files);
    files(f).blockWn = {files(f).prewn files(f).postwn};
    files(f).blockDark = {files(f).predark files(f).postdark};
    files(f).blockBar = {files(f).prebars files(f).postbars};
    files(f).blockDrift = {files(f).predrift files(f).postdrift};
end

