pathname = '\\angie\Angie_analysis\DOI_experiments';
n=0;

% n=n+1;
% files(n).expt = '102815';
% files(n).dir = '10_28_15';
% files(n).tank = '';
% files(n).clusterfile = 'cluster_data_10_28_15_HTR2A_doi_cluster_102815';
% files(n).analysisfile = 'analysis_102815b';
% files(n).predark = 'darrk_predoi1';
% files(n).postdark = 'dark_postdoi1';
% files(n).prewn = 'wn_predoi1';
% files(n).postwn = 'wn_postdoi1';
% files(n).predrift = 'drift_predoi1';
% files(n).postdrift = 'drift_postdoi1';
% files(n).prebars = '';
% files(n).postbars = '';
% files(n).prebackground = 'background_predoi1';
% files(n).postbackground = 'background_postdoi1';
% files(n).prepinpFile = 'laser_int4_5'
% files(n).postpinpFile = 'laser_int4_5_post'
% files(n).injection = 'injection';
% files(n).treatment = 'DOI';
% files(n).layers = [];
% files(n).tip1 = 525; files(n).tip2 = 500; files(n).angle = 45; %not sure need to check
% files(n).notes = 'good data';%'no UDP';
% files(n).misc = 'HTR2A_Pinp';

%%lost units after pinping, put pre sessions to keep compile going
n=n+1;
files(n).expt = '042016';
files(n).dir = '04_20_16';
files(n).tank = '';
files(n).clusterfile = 'cluster_data_04_20_16_cluster_042016';
files(n).analysisfile = 'analysis_042016b';
files(n).predark = 'dark_pre1';
files(n).postdark = 'dark_pre1';
files(n).prewn = 'wn_pre1';
files(n).postwn = 'wn_pre1';  
files(n).predrift = 'drift_pre1';
files(n).postdrift = 'drift_pre1';
files(n).prebars = 'bars_pre1';
files(n).postbars = '';
files(n).injection = '';
files(n).treatment = 'DOI';
files(n).prepinpFile = 'laser_pre_int_9_5'
files(n).postpinpFile = 'laser_pre_int_9_5' %just a test...not a post session
files(n).layers = [];
files(n).notes = 'good data';
files(n).misc = 'HTR2A_Pinp'% few units - lost from laser 
files(n).tip1 = 550; files(n).tip2 = 500; files(n).angle = 45;  
%files(n).wn_camera = {'eyefile1.mat','eyefile2.mat'};  

n=n+1;
files(n).expt = '040716';
files(n).dir = '04_07_16';
files(n).tank = '';
files(n).clusterfile = 'cluster_data_04_07_16_cluster_040716';
files(n).analysisfile = 'analysis_040716b';
files(n).predark = 'dark_pre1';
files(n).postdark = 'dark_post1';
files(n).prewn = 'wn_pre1';
files(n).postwn = 'wn_post1';  
files(n).predrift = 'drift_pre1';
files(n).postdrift = 'drift_post1';
files(n).prebars = 'bars_pre1';
files(n).postbars = 'bars_post3';
files(n).injection = 'injection';
files(n).treatment = 'DOI';
files(n).prepinpFile = 'laser_post_3_5' % just a test
files(n).postpinpFile = 'laser_post_3_5' %pinped post only 
files(n).layers = [];
files(n).notes = 'good data'; 
files(n).misc = 'HTR2A_Pinp';
files(n).tip1 = 525; files(n).tip2 = 500; files(n).angle = 45;  
files(n).predark_camera = 'bars_pre_040716_eye';  
files(n).postdark_camera = 'bars_post_040716_eye';% 
files(n).badsites = [];

n=n+1;
files(n).expt = '010816';
files(n).dir = '01_08_16';
files(n).tank = '';
files(n).clusterfile = 'cluster_data_01_08_16htr2a_doi_cluster_010816';
files(n).analysisfile = 'analysis_010816b';
files(n).predark = 'dark_predoi1';
files(n).postdark = 'dark_postdoi1';
files(n).prewn = 'wn_predoi1';
files(n).postwn = 'wn_postdoi1';
files(n).predrift = 'drift_predoi1';
files(n).postdrift = 'drift_postdoi1';
files(n).prebars = '';
files(n).postbars = '';
files(n).prepinpFile = 'laser_pre_int_075' %psth_power1
files(n).postpinpFile ='laser_postdoi_int075' %psth_power2
files(n).injection = '';
files(n).treatment = 'DOI';
files(n).layers = [];
files(n).tip1 = 525; files(n).tip2 = 500; files(n).angle = 45; %not sure, need to check
files(n).notes = 'good data'; % few units, most drifted after pre pinping
files(n).misc = 'HTR2A_Pinp';

for f = 1:length(files);
    files(f).blockWn = {files(f).prewn files(f).postwn};
    files(f).blockDark = {files(f).predark files(f).postdark};
    files(f).blockBar = {files(f).prebars files(f).postbars};
    files(f).blockDrift = {files(f).predrift files(f).postdrift};
end

% 
% plot(psth_power2(5,:));hold on; plot(psth_power4(5,:))
% plot([50 50],[0 50])