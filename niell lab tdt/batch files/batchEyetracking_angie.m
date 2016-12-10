pathname = '\\angie\Angie_analysis\DOI_experiments';
n=0;

n=n+1;
files(n).expt = '120816';
files(n).dir = '12_08_16_eye';
files(n).Tank_Name = '12_08_16_eye';
files(n).prehigh = 'drift_highcontrast_pre';
files(n).posthigh = 'drift_highcontrast_post_eye';
files(n).prelow = 'drift_lowcontrast_pre';
files(n).postlow = 'drift_lowcontrast_post_eye';  
files(n).prehigh_camera = '12_08_16_drift_highcontrast_pre_eye';  
files(n).posthigh_camera = '12_08_16_drift_highcontrast_post_eye';
files(n).prelow_camera = '12_08_16_drift_lowcontrast_pre_eye';  
files(n).postlow_camera = '12_08_16_drift_lowcontrast_post_eye';
files(n).treatment = 'DOI'; 
files(n).notes = 'good data'; 
files(n).misc = ''; 


for f = 1:length(files);
    files(f).blockHigh = {files(f).prehigh files(f).posthigh};
    files(f).blockLow = {files(f).prelow files(f).postlow};
   % files(f).blockDark = {files(f).predark files(f).postdark};
end
