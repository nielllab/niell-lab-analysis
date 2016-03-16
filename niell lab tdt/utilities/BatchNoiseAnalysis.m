dbstop if error

afile = {   'D:\Jen_analysis\NR5A_Pinping\6_17_15\analysis.mat'};                       

clusterfilenames= { 'D:\Jen_analysis\NR5A_Pinping\6_17_15\cluster_data_6_17_15_6_17_15.mat'};               % MLRstim #16}



for i = 1:length(afile);
    i
    afile
    
    analysisFile = afile{i};
    load(analysisFile);
   
    clusterFile = clusterfilenames{i};
    load(clusterFile);
    
    
    for i = 1:length(Block_Name)
        if strcmp('wn',Block_Name{i}(1:2)) %'sp' if spot analysis 'wn' if cm_noise
            blocknum=i;
        end
    end
    pdfname = [analysisFile(1:end-4) '.pdf'];
    
    
movieFile = 'C:\Users\lab\Desktop\movieFiles\cortex\wn_cortex_012alpha1_5hzLg30Hz.mat';
load(movieFile);
    
noise_analysis(clusterFile,analysisFile,pdfname,movieFile,Block_Name{blocknum},blocknum,1, 1)
% movietype==fl_noise==2; 
% movietype==cm_noise==1 
% movietype==mv_noise==3
 
    end
    