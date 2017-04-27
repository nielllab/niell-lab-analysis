pathname = '\Tai_analysis';
n=0;

n=n+1;
files(n).expt = '041817';
files(n).dir = '04_18_17';
files(n).tank = '';
files(n).clusterfile = 'cluster_data_4_18_17_cluster_041817';
files(n).analysisfile = 'analysis_041817';
files(n).predark = 'dark_pre1';
files(n).postdark = 'dark_post1';
files(n).prewn = '';
files(n).postwn = '';  
files(n).predrift = 'drift_pre1';
files(n).postdrift = 'drift_post1';
files(n).injection = 'injection';
files(n).treatment = 'DOI'; %1 mg/g
files(n).notes = 'good data'; 
files(n).misc = ''; 
files(n).tip1 = 700; files(n).tip2 = 700; files(n).angle = 90;%vertical penetration 
files(n).badsites = [];


n=n+1;
files(n).expt = '041417';
files(n).dir = '04_14_17';
files(n).tank = '';
files(n).clusterfile = 'cluster_data_04_14_17_cluster_041417';
files(n).analysisfile = 'analysis_041417';
files(n).predark = 'dark_pre1';
files(n).postdark = 'dark_post1';
files(n).prewn = 'wn_pre1';
files(n).postwn = 'wn_post1';  
files(n).predrift = '';
files(n).postdrift = '';
files(n).injection = 'injection';
files(n).treatment = 'DOI'; %1 mg/g
files(n).notes = 'good data'; 
files(n).misc = ''; 
files(n).tip1 = 700; files(n).tip2 = 700; files(n).angle = 90;%vertical penetration 
files(n).badsites = [];


n=n+1;
files(n).expt = '041117';
files(n).dir = '04_11_17b';
files(n).tank = '';
files(n).clusterfile = 'cluster_data_04_11_17_cluster_041117';
files(n).analysisfile = 'analysis_041117';
files(n).predark = 'dark_pre1';
files(n).postdark = 'dark_post1';
files(n).prewn = '';
files(n).postwn = '';  
files(n).predrift = 'drift_pre1';
files(n).postdrift = 'drift_post1';
files(n).injection = 'injection';
files(n).treatment = 'Saline'; 
files(n).notes = 'good data'; 
files(n).misc = ''; 
files(n).tip1 = 700; files(n).tip2 = 700; files(n).angle = 90;%vertical penetration 
files(n).badsites = [];



for f = 1:length(files);
    files(f).blockWn = {files(f).prewn files(f).postwn};
    files(f).blockDark = {files(f).predark files(f).postdark};
   % files(f).blockBar = {files(f).prebars files(f).postbars};
    files(f).blockDrift = {files(f).predrift files(f).postdrift};
end