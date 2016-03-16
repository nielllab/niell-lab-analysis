function wn_mv = getWn_mv(clustfile,afile, block,redo,periodFrames)
%%% read in single unit spike times for a given block
%%% this is mostly just a matter of sorting the spikes from one block
%%% but nice to do it just in one line!


load(clustfile,'Block_Name','Tank_Name','frameEpocs');
blocknum = find(strcmp(block,Block_Name));
Block_Name = Block_Name{blocknum}

load(afile,'wn_mv');

if ~exist('wn_mv','var') | length(wn_mv)<blocknum  | isempty(wn_mv(blocknum).frms) | redo
    try
        display('getting spikes')
        spikes=getSpikes(clustfile,afile,block,redo);
        display('getting speed')
        spd = getSpeed(clustfile,afile,block,redo);
        frms = frameEpocs{blocknum}(1,:);
        cycFrames = mod(frms,periodFrames)/periodFrames;
        frameT = frameEpocs{blocknum}(2,:);
        display('doing hist')
        frameSpd = interp1(spd.t,spd.v,frameT,'nearest');
        for c =1 :length(spikes.sp);
            c
            R = histc(spikes.sp{c},frameT);
            R = R(1:end-1)./diff(frameT);
            wn_mv(blocknum).R(c,:) = R;
            for i = 1:20;
                for mv = 1:2;
                    if mv ==1
                        trials = find(frameSpd<1 & cycFrames>(i-1)/20 & cycFrames<=i/20);
                    else
                        trials = find(frameSpd>1 & cycFrames>(i-1)/20 & cycFrames<=i/20);
                    end
                    trials = trials(trials<=length(R));
                    crf(i,mv) = mean(R(trials));
                    crf_err(i,mv) = std(R(trials))/sqrt(length(trials));
                end
            end
             wn_mv(blocknum).crf(c,:,:) = crf;
            wn_mv(blocknum).crf_err(c,:,:) = crf_err;
            wn_mv(blocknum).spont(c,:)= mean(crf([1 20],:),1);
            wn_mv(blocknum).evoked(c,:) = mean(crf(9:12,:),1) - wn_mv(blocknum).spont(c,:);
            
            [r m b] = regression(crf(:,1)'- mean(crf([1 20],1),1),crf(:,2)'- mean(crf([1 20],1),1));
            wn_mv(blocknum).gain(c) =m;
            wn_mv(blocknum).intercept(c) = b;
            wn_mv(blocknum).reliability(c)=r;
            figure
            subplot(1,2,1)
            errorbar(crf,crf_err); xlim([0.5 20.5])
            subplot(1,2,2)
            plot(crf(:,1)- mean(crf([1 20],1),1),crf(:,2)- mean(crf([1 20],1),1),'o')  ; hold on; plot([0 5],[0 5]);
            plot([0 5],[0 5]*m +b,'b');

        end
        wn_mv(blocknum).frms = frms(1:end-1);
        wn_mv(blocknum).cycFrames = cycFrames(1:end-1);
        wn_mv(blocknum).frameSpd = frameSpd(1:end-1);
        
    catch
        wn_mv(blocknum).frms = [];
        wn_mv(blocknum).cycFrames = [];
        wn_mv(blocknum).frameSpd = [];
        wn_mv(blocknum).crf =[];
        wn_mv(blocknum).R = [];
    end
    save(afile,'wn_mv','-append')
end

wn_mv = wn_mv(blocknum);
