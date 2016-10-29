function wn_mv = getWn_mv(clustfile,afile, block,redo,periodFrames)
%%% read in single unit spike times for a given block
%%% this is mostly just a matter of sorting the spikes from one block
%%% but nice to do it just in one line!


load(clustfile,'Block_Name','Tank_Name','frameEpocs');
blocknum = find(strcmp(block,Block_Name));
Block_Name = Block_Name{blocknum}

load(afile,'wn_mv');
load(afile,'wn');

if ~exist('wn_mv','var') | length(wn_mv)<blocknum  | ~isfield(wn_mv(blocknum),'frameR')  | isempty(wn_mv(blocknum).frameR) | redo
   % try
        display('getting spikes')
        spikes=getSpikes(clustfile,afile,block,0);
        display('getting speed')
        spd = getSpeed(clustfile,afile,block,0);
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
            for f = 1:max(frms);
                for mv=1:2
                    if mv==1
                        trials = find(frameSpd<0.5 & frms ==f);
                    else
                        trials = find(frameSpd>0.5 & frms ==f);
                    end
                    trials = trials(trials<=length(R));
                    frameR(f,mv) = mean(R(trials));
                end
            end
            wn_mv(blocknum).frameR(c,:,:) = frameR;
            
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
%             figure
%             subplot(1,2,1)
%             errorbar(crf,crf_err); xlim([0.5 20.5])
%             subplot(1,2,2)
%             plot(crf(:,1)- mean(crf([1 20],1),1),crf(:,2)- mean(crf([1 20],1),1),'o')  ; hold on; plot([0 5],[0 5]);
%             plot([0 5],[0 5]*m +b,'b');

        end
       
%          if exist('params','var');
%             
%             
%             all_img_STA(cellrange)= all_img;
%             all_fit_STA(cellrange)=all_fit;
%             STA_nx(cellrange)=field2array(params,'nx');
%             STA_ny(cellrange)=field2array(params,'ny');
%             STA_phase(cellrange)=field2array(params,'phase');
%             STA_sigx(cellrange)=field2array(params,'sigx');
%             STA_sigy(cellrange)=field2array(params,'sigy');
%             STA_exp_var(cellrange)=field2array(params,'exp_var');
%             
%             
%         else
%             STA_nx(cellrange)= NaN;
%             STA_ny(cellrange)=  NaN;
%             STA_phase(cellrange)= NaN;
%             STA_sigx(cellrange)= NaN;
%             STA_sigy(cellrange)= NaN;
%             STA_exp_var(cellrange)=NaN;
%             
%          end
%          
%          clear m ind x y t_lag STA STA1 STA_post STA1_post c crf crf1 crf_post crf1_post
%          
%          if exist ('wn','var')
%              for w = 1:length(wn)
%                  sta[s,t,E] = sta(spikes.sp,lfpMove,smp,plt,w,T,D,err)
%                 wn(blocknum).STA = sta;
%                  
%                  %%%Dtermine time point with maximial response
%                  [m ind] = max(abs(STA(:)-127));
%                  [x y t_lag] = ind2sub(size(STA),ind);
%                  
%                  STA1{w} = STA(:,:,t_lag)-128;
%                  figure
%                  imagesc(STA1{1,w}',[-64 64]); axis equal
%              end
%          end
%              STA_peak(cellrange)=STA1
            
        wn_mv(blocknum).frms = frms(1:end-1);
        wn_mv(blocknum).cycFrames = cycFrames(1:end-1);
        wn_mv(blocknum).frameSpd = frameSpd(1:end-1);
        
%     catch
%         wn_mv(blocknum).frms = [];
%         wn_mv(blocknum).cycFrames = [];
%         wn_mv(blocknum).frameSpd = [];
%         wn_mv(blocknum).crf =[];
%         wn_mv(blocknum).R = [];
%     end
    save(afile,'wn_mv','-append')
end

wn_mv = wn_mv(blocknum);
