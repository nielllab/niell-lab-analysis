clear all

n=0;


afile = {    'C:\data\lgn rf project\lgn_analysis\030912_rec3_analysis3_final.mat',...
    'C:\data\lgn rf project\lgn_analysis\030912_rec4_analysis_4final.mat',...
    'C:\data\lgn rf project\lgn_analysis\031012_analysis1check.mat',...
    'C:\data\lgn rf project\lgn_analysis\031012_analysis2check.mat',...
    'C:\data\lgn rf project\lgn_analysis\031512_analysis2check.mat',...
    'C:\data\lgn rf project\lgn_analysis\031512_analysis3check.mat',...
    'C:\data\lgn rf project\lgn_analysis\031612_analysis1checkb.mat',...
    'C:\data\lgn rf project\lgn_analysis\031612_analysischeck2b.mat',...
    'C:\data\lgn rf project\lgn_analysis\031612_analysis3check.mat', ...
    'C:\data\lgn rf project\lgn_analysis\052512_analysis1check.mat',...
    'C:\data\lgn rf project\lgn_analysis\060112_analysis1check.mat',...
    'C:\data\lgn rf project\lgn_analysis\060712_analysis1check.mat',...
    'C:\data\lgn rf project\lgn_analysis\060812_analysis1check.mat',...
    'C:\data\lgn rf project\lgn_analysis\061212_analysis1check.mat',...
    'C:\data\lgn rf project\lgn_analysis\061212_analysis2check.mat',...
    'C:\data\lgn rf project\lgn_analysis\062212_analysis3check.mat',...
    'C:\data\lgn rf project\lgn_analysis\062212_analysis4check.mat',...
    'C:\data\lgn rf project\lgn_analysis\071612_analysis1new.mat', ...
    'C:\data\lgn rf project\lgn_analysis\071712_analysis1new.mat' ...
    'C:\data\lgn rf project\lgn_analysis\071712_analysis2new.mat' ...
    'C:\data\lgn rf project\lgn_analysis\071812_analysis1new.mat' ...
    'C:\data\lgn rf project\lgn_analysis\071812_analysis2new.mat' ...
    'C:\data\lgn rf project\lgn_analysis\072012_analysis1newst.mat' ...
    'C:\data\lgn rf project\lgn_analysis\072012_analysis2new.mat' ...
    'C:\data\lgn rf project\lgn_analysis\072012_analysis3.mat' ...
    'C:\data\lgn rf project\lgn_analysis\072412_analysis1.mat' ...
    'C:\data\lgn rf project\lgn_analysis\072512_analysis1.mat' ...
    'C:\data\lgn rf project\lgn_analysis\072512_ANALYSIS2.mat' ...
    'C:\data\lgn rf project\lgn_analysis\080212_analysis1.mat' ...
    'C:\data\lgn rf project\lgn_analysis\080212_analysis3.mat' ...
    'C:\data\lgn rf project\lgn_analysis\080312_analysis2.mat' ...
    'C:\data\lgn rf project\lgn_analysis\080312_analysis3bnew.mat' };



for fnum = 1:length(afile)
    fnum
    analysisFile = afile{fnum}
    load(analysisFile);
    clusterfilename
    pname
    if exist(clusterfilename,'file')
        clusterFile = clusterfilename;
    elseif exist(clusterfilename((length(pname)+1):end),'file')
        clusterFile = clusterfilename((length(pname)+1):end);
    elseif exist([clusterfilename((length(pname)+1):end) '.mat'],'file')
        clusterFile = [clusterfilename((length(pname)+1):end) '.mat'];
    else
        [fname pname] = uigetfile('*.mat','cluster file');
        clusterFile = fullfile(pname,fname);
        clusterfilename = clusterFile;
        if fname~=0
            save(analysisFile,'clusterfilename','-append');
        end
    end
    clusterFile
    load(clusterFile);
    blocknum=0;
    for b = 1:length(Block_Name)
        if strcmpi('wn',Block_Name{b}(1:2))
            blocknum=b;
        end
    end
    
    blocknum
    Block_Name{blocknum}
    flags = struct('visStim',1);
    
    
    data = getTDTdata(Tank_Name,Block_Name{blocknum},1:4:max(cells(:,1)),flags);
    if isnan(data.frameEpocs)
        sprintf('check if tank %s is registered',Tank_Name)
        break
    end
    % blockdata.mouseT = data.mouseT;
    % blockdata.mouseV = data.mouseV;
    exptdata(fnum).block = Block_Name{blocknum};
    exptdata(fnum).tank = Tank_Name;
    exptdata(fnum).analysis_file = afile;
    exptdata(fnum).cluster_file = clusterFile;
    exptdata(fnum).stimEpocs=data.stimEpocs;
    exptdata(fnum).frameEpocs = data.frameEpocs;
    
    for cell_n = 1:size(cells,1);
        unitdata(cell_n+n).spikes =spikeT{cell_n}(spikeT{cell_n}>(blocknum-1)*10^5 & spikeT{cell_n}<(blocknum-1 + 0.5)*10^5) - (blocknum-1)*10^5;
        unitdata(cell_n+n).expnum = fnum;
        %     blockdata.lfpT{cell_n} = data.lfpT{cells(cell_n,1)};
        %     blockdata.lfpV{cell_n} = data.lfpData{cells(cell_n,1)};
    end
    n=n+size(cells,1);
end

manual_type = xlsread('C:\data\lgn rf project_new0824\lgn_analysis\lgn types 091412.xlsx','A1:A294');
manual_outside = xlsread('C:\data\lgn rf project_new0824\lgn_analysis\lgn types 091412.xlsx','B1:B294');
manual_outside(isnan(manual_outside))=0;

used =(manual_type~=0);
unitdata = unitdata(used);


load('C:\data\movie files\wn016alpha1_10hzLg60Hz.mat');

save('Piscopo_wn_data072513','exptdata','unitdata','moviedata');