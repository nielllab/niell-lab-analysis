function speed = getSpeed(clustfile,afile, block,redo )
%%% read in optical mouse speed data for a given block
%%% checks to see if data already exists
%%% otherwise reads it in and saves to analysis file, with entry for each block
%%% returns data for specified block
%%% speed.t = timestamps, speed.v = speed


load(clustfile,'Block_Name','Tank_Name');
blocknum = find(strcmp(block,Block_Name));
Block_Name = Block_Name{blocknum}

load(afile,'speed');

if ~exist('speed','var') | length(speed)<blocknum  | isempty(speed(blocknum).t) | redo
    try
        flags = struct('mouseOn',1);
        tdtData= getTDTdata(Tank_Name, Block_Name, 1, flags);
        speed(blocknum).t= tdtData.mouseT;
        speed(blocknum).v = tdtData.mouseV;
    catch
        speed(blocknum).t=[];
        speed(blocknum).v=[];
    end
    save(afile,'speed','-append')
end

speed = speed(blocknum);
