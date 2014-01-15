afile = {...
%    'C:\data\ephys matlab data\021412_awake_chr2\wndrift1\analysis.mat', ...     %%LFP problems
%     'C:\data\ephys matlab data\021512_awake_pptg\wn1b\analysis.mat', ...        %%LFP problems
    'C:\data\ephys matlab data\021612_awake\wn1drift1\analysis.mat', ...          % MLRstim #1; great desynchronization
    'C:\data\ephys matlab data\021612_awake\wn2b\analysis.mat',...                % MLRstim #2; good desychronization
    'C:\data\ephys matlab data\021612_awake\wn3\analysis.mat',...                 % MLRstim #3; good desychronization at 50Hz
    'C:\data\ephys matlab data\021612_awake\wn4bdrift5\analysis.mat',...          % MLRstim #4    
    'C:\data\ephys matlab data\021812_awake_pptg\wn2\analysis.mat',...            % MLRstim #5    
    'C:\data\ephys matlab data\021812_awake_pptg\wn3e_drift2\analysis.mat',...    % MLRstim #6
    'C:\data\ephys matlab data\021812_awake_pptg\wn5\analysis.mat',...            % MLRstim #7
    'C:\data\ephys matlab data\070112_awake_mlr\wn1fg\analysis.mat',...           % MLRstim #8
    'C:\data\ephys matlab data\070112_awake_mlr\wn2a\analysis.mat',...            % MLRstim #9
    'C:\data\ephys matlab data\070112_awake_mlr\wn3d\analysis.mat',...            % MLRstim #10
    'C:\data\ephys matlab data\070112_awake_mlr\wn4\analysis.mat',...             % MLRstim #11  
    'C:\data\ephys matlab data\070112_awake_mlr\wn5b\analysis.mat',...            % MLRstim #12
    'C:\data\ephys matlab data\070312_awake_mlr\wn1h\analysis.mat',...            % MLRstim #13
    'C:\data\ephys matlab data\070312_awake_mlr\wn2e\analysis.mat',...            % MLRstim #14
    'C:\data\ephys matlab data\070312_awake_mlr\wn3a\analysis.mat',...            % MLRstim #15
    'C:\data\ephys matlab data\070312_awake_mlr\wn4c\analysis.mat'}               % MLRstim #16
tic

for i = 14:length(afile);
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
    
    