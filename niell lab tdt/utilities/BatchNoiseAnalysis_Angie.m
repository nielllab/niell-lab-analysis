dbstop if error

afile = {'D:\Angie_analysis\DOI_experiments\03_09_16\analysis_030916.mat'};                       

clusterfilenames= { 'D:\Angie_analysis\DOI_experiments\03_09_16\cluster_data_03_09_16_cluster_030916.mat'}
             % MLRstim #16}



for i = 1:length(afile);
    i
    afile
    
    analysisFile = afile{i};
    load(analysisFile);
   
    clusterFile = clusterfilenames{i};
    load(clusterFile);
    
%     [fname pname] = uigetfile('*.mat','cluster file');
%      clusterFile = fullfile(pname,fname);
%      clusterfilename = clusterFile;     
%      save(analysisFile,'clusterfilename','-append');
  
  
   
    for i = 1:length(Block_Name)
        if strcmp('sp',Block_Name{i}(1:2)) %'sp' if spot analysis 'wn' if cm_noise
            blocknum=i;
        end
    end
    pdfname = [analysisFile(1:end-4) '.pdf'];
    
    
movieFile = 'C:\Users\Angie Michaiel\Desktop\movie files\cortex\wn_cortex_012alpha1_5hzLg30Hz.mat';
load(movieFile);
    

    %pptgSpeedAnalysis(clusterFile,analysisFile,pdfname,Block_Name{blocknum},blocknum,1, 30);
    noise_analysis(clusterFile,analysisFile,pdfname,movieFile,Block_Name{blocknum},blocknum,2, 1)
% movietype==fl_noise==2; 
% movietype==cm_noise==1 
% movietype==mv_noise==3
    toc

    %plot(vdata)
    end
    