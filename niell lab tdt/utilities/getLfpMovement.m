function lfpMove = getLfpMovement(clustfile,afile,block, redo)
%%% read in lfp spectrum data for a given block
%%% checks to see if data already exists
%%% otherwise reads it in and saves to analysis file, with entry for each block
%%% returns data for specified block
%%% lfp.t = timestamps, lfp.freq = frequencies 
%%% lfp.normspect = spectrum for each channel


load(clustfile,'Block_Name','Tank_Name');
blocknum = find(strcmp(block,Block_Name));
Block_Name = Block_Name{blocknum}

load(afile,'lfpMove');

if ~exist('lfpMove','var') | length(lfpMove)<blocknum  | isempty(lfpMove(blocknum).freq) | redo
  
   try
   lfp = getLFP(clustfile,afile,block,0);
   spd = getSpeed(clustfile,afile,block,0);
   s = lfp.normspect;
   v = interp1(spd.t,spd.v,lfp.t);
   for ch = 1:size(s,1);
       for i = 1:2
           if i==1
               meanS =  squeeze(mean(s(ch,v<=1,:),2));
           else
               meanS =  squeeze(mean(s(ch,v>1,:),2));
           end
       meanSpect(ch,:,i) = interp1(lfp.freq,meanS,0.5:0.5:100);
       end
   end
   lfpMove(blocknum).meanSpect = meanSpect;
   lfpMove(blocknum).freq = 0.5:0.5:100;
   %keyboard
   catch ME
        lfpMove(blocknum).meanSpect=[];
        lfpMove(blocknum).freq=[];
        display('couldnt get lfp movement')
       % getReport(ME)
    end
    save(afile,'lfpMove','-append')
end

lfpMove = lfpMove(blocknum);





