function lfp = getLFP(clustfile,afile,block, redo)
%%% read in lfp spectrum data for a given block
%%% checks to see if data already exists
%%% otherwise reads it in and saves to analysis file, with entry for each block
%%% returns data for specified block
%%% lfp.t = timestamps, lfp.freq = frequencies 
%%% lfp.normspect = spectrum for each channel


load(clustfile,'Block_Name','Tank_Name','ch');
nChan=ch;
blocknum = find(strcmp(block,Block_Name));
Block_Name = Block_Name{blocknum}

load(afile,'lfp');


if ~exist('lfp','var') | length(lfp)<blocknum  | isempty(lfp(blocknum).t) | redo
   try
        flags = struct('lfpTseries',1,'lfpSpectra',1);
        tdtData= getTDTdata(Tank_Name, Block_Name, 1:nChan, flags);
        lfp(blocknum).freq = tdtData.spectF{1};
        lfp(blocknum).t = tdtData.spectT{1};
        for ch = 1:nChan;
            lfpchan = tdtData.spectData{ch};
            normalizer = 1:size(lfpchan,2);
            normalizer = repmat(normalizer,size(lfpchan,1),1);
            normspect(ch,:,:) = lfpchan.*normalizer;
        end
        %%% downsample lfp
        for f = 1:nChan/4;
            ns(f,:,:) = median(normspect((f-1)*4 + (1:4),:,4:4:end),1);
        end
        lfp(blocknum).normspect = ns;
        lfp(blocknum).freq = lfp(blocknum).freq(4:4:end);
    catch
        lfp(blocknum).t=[];
        lfp(blocknum).freq=[];
        lfp(blocknum).normspect=[];
    end

    save(afile,'lfp','-append')
end

lfp = lfp(blocknum);





