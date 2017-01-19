pathname = '\\angie\Angie_analysis\DOI_experiments';
n=0;


n=n+1;
files(n).expt = '011817';
files(n).dir = '01_18_17_eye';
files(n).Tank_Name = '01_18_17_eye';
files(n).predetection = 'detection_pre1';
files(n).postdetection = 'detection_post1';
files(n).predetection_camera='01_18_17_detection_pre_eye';
files(n).postdetection_camera='01_18_17_detection_post_eye';
files(n).treatment = 'DOI'; 
files(n).notes = 'good data'; 
files(n).misc = 'J307c'


n=n+1;
files(n).expt = '011317';
files(n).dir = '01_13_17_eye';
files(n).Tank_Name = '01_13_17_eye';
files(n).predetection = 'detection_pre1';
files(n).postdetection = 'detection_post2';
files(n).predetection_camera='01_13_17_detection_pre_eye';
files(n).postdetection_camera='01_13_17_detection_post_eye';
files(n).treatment = 'Saline'; 
files(n).notes = 'good data'; 
files(n).misc = 'J294c'

n=n+1;
files(n).expt = '011117';
files(n).dir = '01_11_17_eye';
files(n).Tank_Name = '01_11_17';
files(n).predetection = 'detection-pre1';
files(n).postdetection = 'detection_post2';
files(n).predetection_camera='01_11_17_detection_pre_eye';
files(n).postdetection_camera='01_11_17_detection_post_eye';
files(n).treatment = 'HT'; 
files(n).notes = 'good data'; 
files(n).misc = 'J294e'%with ephys


% n=n+1;
% files(n).expt = '010417';
% files(n).dir = '01_04_17_eye';
% files(n).Tank_Name = '01_04_17';
% files(n).predetection = 'detection_pre1';
% files(n).postdetection = 'detection_post1';
% files(n).predetection_camera='01_04_17_detection_pre_eye';
% files(n).postdetection_camera='01_04_17_detection_post_eye';
% files(n).treatment = 'Saline'; 
% files(n).notes = 'good data'; 
% files(n).misc = 'J294c'

n=n+1;
files(n).expt = '122216';
files(n).dir = '12_22_16_eye';
files(n).Tank_Name = '12_22_16_eye';
files(n).predetection = 'detection_pre1_eye';
files(n).postdetection = 'detection_post_eye';
files(n).predetection_camera='12_22_16_detection_pre_eye';
files(n).postdetection_camera='12_22_16_detection_post_eye';
files(n).treatment = 'Saline'; 
files(n).notes = 'good data'; 
files(n).misc = 'J294b'


n=n+1;
files(n).expt = '121916';
files(n).dir = '12_19_16a_eye';
files(n).Tank_Name = '12_19_16a_eye';
files(n).predetection = 'detection_pre2_eye';
files(n).postdetection = 'detection_post1';
files(n).predetection_camera='12_19_16a_detection_pre_eye';
files(n).postdetection_camera='12_19_16a_detection_post_eye';
files(n).treatment = 'DOI'; 
files(n).notes = 'good data'; 
files(n).misc = 'J294b'; %second time with DOI (11 days in between)
                         % -- behaviorally didn't look like DOI was as strong as first time

    %stim started after cameraTTL (first time starting this way)
               

n=n+1;
files(n).expt = '121316';
files(n).dir = '12_13_16_eye';
files(n).Tank_Name = '12_13_16_eye';
files(n).predetection = 'detection_pre1_eye';
files(n).postdetection = 'detection_post_eye';
files(n).predetection_camera='12_13_16_detection_pre_eye';
files(n).postdetection_camera='12_13_16_detection_post_eye';
files(n).treatment = 'DOI'; 
files(n).notes = 'good data'; 
files(n).misc = 'J294c'; %first time with DOI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OLD STIM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n=n+1;
% files(n).expt = '121216';
% files(n).dir = '12_12_16_eye';
% files(n).Tank_Name = '12_12_16_eye';
% files(n).prehigh = 'drift_highcontrast3_pre_eye';
% files(n).posthigh = 'drift_highcontrast_post_eye2';
% files(n).prelow = 'drift_lowcontrast2_pre_eye';
% files(n).postlow = 'drift_lowcontrast_post_eye';  
% files(n).prehigh_camera = '12_12_16_drift_highcontrast_pre_eye';  
% files(n).posthigh_camera = '12_12_16_drift_highcontrast_post_eye';
% files(n).prelow_camera = '12_12_16_drift_lowcontrast_pre_eye';  
% files(n).postlow_camera = '12_12_16_drift_lowcontrast_post_eye';
% files(n).treatment = 'Saline'; 
% files(n).notes = 'good data'; 
% files(n).misc = 'J294c'; 
% 
% 
% n=n+1;
% files(n).expt = '120816';
% files(n).dir = '12_08_16_eye';
% files(n).Tank_Name = '12_08_16_eye';
% files(n).prehigh = 'drift_highcontrast_pre';
% files(n).posthigh = 'drift_highcontrast_post_eye';
% files(n).prelow = 'drift_lowcontrast_pre';
% files(n).postlow = 'drift_lowcontrast_post_eye';  
% files(n).prehigh_camera = '12_08_16_drift_highcontrast_pre_eye';  
% files(n).posthigh_camera = '12_08_16_drift_highcontrast_post_eye';
% files(n).prelow_camera = '12_08_16_drift_lowcontrast_pre_eye';  
% files(n).postlow_camera = '12_08_16_drift_lowcontrast_post_eye';
% files(n).treatment = 'DOI'; 
% files(n).notes = 'good data'; 
% files(n).misc = 'J294b'; 


for f = 1:length(files);
    files(f).blockDetect = {files(f).predetection files(f).postdetection};
    %     files(f).blockHigh = {files(f).prehigh files(f).posthigh};
    %     files(f).blockLow = {files(f).prelow files(f).postlow};
    %     files(f).blockDark = {files(f).predark files(f).postdark};
end
