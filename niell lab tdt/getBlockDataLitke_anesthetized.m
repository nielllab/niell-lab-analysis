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
    %for fnum=1:2
    fnum
    analysisFile = afile{fnum};
    load(analysisFile);
    clusterfilename;
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
    for blocknum = 1:length(Block_Name)
        
        
        blocknum
        Block_Name{blocknum};
        
        flags = struct('visStim',1,'lfpTseries',1,'lfpSpectra',1);
        % flags = struct('visStim',1);
        
        ch = round(median(cells(:,1)))
        
        data = getTDTdata(Tank_Name,Block_Name{blocknum},ch,flags);
        if isnan(data.frameEpocs)
            sprintf('check if tank %s is registered',Tank_Name)
            break
        end
        % exptdata{fnum}.mouseT{blocknum} = data.mouseT;
        % blockdata{fnum}.mouseV{blocknum} = data.mouseV;
        exptdata{fnum}.block{blocknum} = Block_Name{blocknum};
        exptdata{fnum}.tank = Tank_Name;
        exptdata{fnum}.analysis_file = afile{fnum};
        exptdata{fnum}.cluster_file = clusterFile;
        exptdata{fnum}.stimEpocs{blocknum}=data.stimEpocs;
        exptdata{fnum}.frameEpocs{blocknum} = data.frameEpocs;
        exptdata{fnum}.lfpT{blocknum} = data.lfpT{ch};
        exptdata{fnum}.lfpV{blocknum} = data.lfpData{ch};
        exptdata{fnum}.displayHeight = displayHeight;
        exptdata{fnum}.displayOffset = displayOffset;
        for cell_n = 1:size(cells,1);
            unitdata{cell_n+n}.spikes{blocknum} =spikeT{cell_n}(spikeT{cell_n}>(blocknum-1)*10^5 & spikeT{cell_n}<(blocknum-1 + 0.5)*10^5) - (blocknum-1)*10^5;
            unitdata{cell_n+n}.expnum = fnum;

        end
    end
    n=n+size(cells,1);
end

% manual_type = xlsread('C:\data\lgn rf project\lgn_analysis\lgn types 091412.xlsx','A1:A294');
% manual_outside = xlsread('C:\data\lgn rf project\lgn_analysis\lgn types 091412.xlsx','B1:B294');
% manual_outside(isnan(manual_outside))=0;
% 
% used =(manual_type~=0);

load 'C:\data\lgn rf project\lgn_analysis\LGNunits'
use_list=find(used);
 for i = 1:length(use_list);
     unitdata_select{i} = unitdata{use_list(i)};
 end
 unitdata=unitdata_select;
 
% 
 save('Piscopo_wn_data090513','exptdata','unitdata');