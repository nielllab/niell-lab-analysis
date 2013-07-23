
[fname, pname] = uigetfile('*.mat','cluster data');
clustfile=fullfile(pname,fname);
load(clustfile);
AllBlocks = Block_Name;

%blocknum = listdlg('ListString',Block_Name,'SelectionMode','single');

[afname, apname] = uigetfile('*.mat','analysis data');
afile = fullfile(apname,afname);
load(afile);
for blocknum=1:length(AllBlocks);
    
    Block_Name = AllBlocks{blocknum}
    
    %flags = struct('visStim',1,'lfpTseries',1,'mouseOn',1,'lfpSpectra',1);
    %flags = struct('visStim',1,'lfpTseries',1,'lfpSpectra',1);
    
    flags = struct('visStim',1);
    
    
    data = getTDTdata(Tank_Name,Block_Name,1:4:max(cells(:,1)),flags);
    
    % blockdata.mouseT = data.mouseT;
    % blockdata.mouseV = data.mouseV;
    blockdata.block = Block_Name;
    blockdata.tank = Tank_Name;
    blockdata.analysis_file = afile;
    blockdata.cluster_file = clustfile;
    blockdata.stimEpocs=data.stimEpocs;
    blockdata.frameEpocs = data.frameEpocs;
    
    for cell_n = 1:size(cells,1);
        blockdata.spikes{cell_n} =spikeT{cell_n}(spikeT{cell_n}>(blocknum-1)*10^5 & spikeT{cell_n}<(blocknum-1 + 0.5)*10^5) - (blocknum-1)*10^5;
        %     blockdata.lfpT{cell_n} = data.lfpT{cells(cell_n,1)};
        %     blockdata.lfpV{cell_n} = data.lfpData{cells(cell_n,1)};
    end
end
