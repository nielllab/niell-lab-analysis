afile = {...
    'C:\data\ephys matlab data\021712_awake\wn1b\analysis.mat', ...         % analyzed
    'C:\data\ephys matlab data\063012_awake\wn3\analysis.mat',...           % reading error
    'C:\data\ephys matlab data\063012_awake\wn5\analysis.mat',...           % analyzed
    'D:\Moses_ephys_data\082712_awake_MLR\wn1a\analysis.mat',...            % analyzed    
    'D:\Moses_ephys_data\082812_awake_MLR\wn1b\analysis.mat',...            % analyzed    
    'D:\Moses_ephys_data\082812_awake_MLR\wn1e\analysis.mat'}               % analyzed
tic

clusterfilenames= {...
    'C:\data\ephys matlab data\021712_awake\wn1b\cluster_data_021712_awake_wn1b.mat', ...          % MLRstim #1; great desynchronization
    'C:\data\ephys matlab data\063012_awake\wn3\cluster_data_063012_awake_mlr_wn3a.mat',...                % MLRstim #2; good desychronization
    'C:\data\ephys matlab data\063012_awake\wn5\cluster_data_063012_awake_mlr_wn5.mat',...                 % MLRstim #3; good desychronization at 50Hz
    'D:\Moses_ephys_data\082712_awake_MLR\wn1a\cluster_data_082712_awake_MLR_wn1a.mat',...          % MLRstim #4    
    'D:\Moses_ephys_data\082812_awake_MLR\wn1b\cluster_data_082812_awake_MLR_wn1b.mat',...            % MLRstim #5    
    'D:\Moses_ephys_data\082812_awake_MLR\wn1e\cluster_data_082812_awake_MLR_wn1e.mat'}               % MLRstim #16}

for i = 3:length(afile);
    i


    analysisFile = afile{i};
    load(analysisFile);

    clusterFiles = clusterfilenames{i};
    load(clusterFiles);
    for i = 1:length(Block_Name)
        if strcmp('wn',Block_Name{i}(1:2))
            blocknum=i;
        end
    end
    pdfname = [analysisFile(1:end-4) '.pdf'];
    
    
    pptgAnalysis(clusterFiles,analysisFile,pdfname,Block_Name{1},1,1, 30);
    toc
    plot(vdata)
    end
    
    