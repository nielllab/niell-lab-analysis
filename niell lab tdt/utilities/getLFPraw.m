function lfpRaw = getLFPraw(clustfile,afile,block, redo)
%%% read in raw lfp traces
%%% checks to see if data already exists
%%% otherwise reads it in and saves to analysis file, with entry for each block
%%% returns data for specified block
%%% lfp.t = timestamps, lfp.freq = frequencies 
%%% lfp.normspect = spectrum for each channel


load(clustfile,'Block_Name','Tank_Name','ch');
nChan=ch;
blocknum = find(strcmp(block,Block_Name));
Block_Name = Block_Name{blocknum}

load(afile,'lfpRaw');


if ~exist('lfpRaw','var') | length(lfpRaw)<blocknum  | isempty(lfpRaw(blocknum).t) | redo
   try
        flags = struct('lfpTseries',1,'lfpSpectra',1)
        tdtData= getTDTdata(Tank_Name, Block_Name, 1:nChan, flags);
        lfpRaw(blocknum).t = tdtData.lfpT{1};
        for f = 1:nChan
            data(:,f) = tdtData.lfpData{f};
        end
        for f = 1:nChan/4;
           d(:,f) = median(data(:,(f-1)*4 + (1:4)),2);
        end
        lfpRaw(blocknum).data = d;
    catch
        lfpRaw(blocknum).t=[];
        lfpRaw(blocknum).data=[];
        
    end

    save(afile,'lfpRaw','-append')
end

lfpRaw = lfpRaw(blocknum);





