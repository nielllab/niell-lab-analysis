afile = {'C:\data\ephys matlab data\070612_awake_mlr\wn4b\analysis.mat',...    %MLR2BFstim#1
    'C:\data\ephys matlab data\070612_awake_mlr\wn5b\analysis.mat',...      %MLR2BFstim#2
    'C:\data\ephys matlab data\070612_awake_mlr\wn6b\analysis.mat',...      %MLR2BFstim#3
    'C:\data\ephys matlab data\070612_awake_mlr\wn7a\analysis.mat',...      %MLR2BFstim#4
    'C:\data\ephys matlab data\070612_awake_mlr\wn8b\analysis.mat',...      %MLR2BFstim#4
    'C:\data\ephys matlab data\070612_awake_mlr\wn9a\analysis.mat',...      %MLR2BFstim#5
    'C:\data\tdt tanks\071312_awake_MLR\wn1b\analysis.mat',...              %MLR2BFstim#6
    'C:\data\ephys matlab data\071312_awake_mlr\wn2d\analysis.mat',...      %MLR2BFstim#7
    'C:\data\ephys matlab data\071312_awake_mlr\wn5\analysis.mat',...       %MLR2BFstim#8
    'D:\Moses_ephys_data\082012_awake_mlr\wn2a\analysis.mat',...            % start here MLR2BFstim#9
    'D:\Moses_ephys_data\082012_awake_mlr\wn3b\analysis.mat',...            %MLR2BFstim#10
    'D:\Moses_ephys_data\082112_awake_mlr\wn2c\analysis.mat',...            %MLR2BFstim#11
    'D:\Moses_ephys_data\082312_awake_MLR\wn3d\analysis.mat',...            %MLR2BFstim#12
    'D:\Moses_ephys_data\082312_awake_MLR\wn5a\analysis.mat',...            %MLR2BFstim#13
    'D:\Moses_ephys_data\082312_awake_MLR\wn6b\analysis.mat'}               %MLR2BFstim#14

for i = 11:length(afile);
    i
close all

    analysisFile = afile{i};
    load(analysisFile);
    clusterFile = clusterfilename((length(pname)+1):end);
    load(clusterFile);
    for i = 1:length(Block_Name)
        if strcmp('wn',Block_Name{i}(1:2))
            blocknum=i;
        end
    end
    pdfname = [analysisFile(1:end-4) '.pdf'];
    
    
    pptgAnalysis(clusterFile,analysisFile,pdfname,Block_Name{blocknum},blocknum,1, 30);
    toc
    end
    
    