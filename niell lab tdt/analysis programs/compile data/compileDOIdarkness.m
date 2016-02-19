clear all
close all

batchDOIephys; %%% load batch file

%%% select the sessions you want based on filters
%%% example
use =  find(strcmp({files.treatment},'DOI') & strcmp({files.notes},'good data') & ~cellfun(@isempty,{files.predark})   )

%use =  find( strcmp({files.treatment},'DOI') & strcmp({files.notes},'good data') & ~cellfun(@isempty,{files.predark}) & ~cellfun(@isempty,{files.postdark}) )

sprintf('%d selected sessions',length(use))

cv2all=[]; %%% initialize data to be saved
meanRall = [];
preCorrs = []; postCorrs=[];
layers = [];

savePDF=0;
redo = 1;
for i = 1:2
    afile = [pathname '\' files(use(i)).dir '\' files(use(i)).analysisfile '.mat'];
    clustfile = [pathname '\' files(use(i)).dir '\' files(use(i)).clusterfile '.mat'] ;
    
    prepostDOIdarkness
    %save(afile,'R','cv2','preCorr','postCorr','-append');
    
    preCorrAll{i}= preCorr; postCorrAll{i}=postCorr;
    preCorrs = [preCorrs preCorr(:)]; postCorrs = [postCorrs postCorr(:)];
    
    
    cv2all = [cv2all; cv2];
    meanRall = [meanRall ; squeeze(mean(R,2))];
    
    eigsFull = zeros(100,2);
    eigsFull(1:size(eigs,1),:)=eigs;
    eigsAll(i,:,:) = eigsFull;
   
    layers = [layers layer];
end



keyboard

  darknessAnalysis
    %%% save out any important data
    prespikes{end+1:end+length(blockSpike)} = blockSpike;
    figure
    plot(blockSpike)
    hold on
    
    blocknm = files(use(i)).postdark;  %%% run for post
    darknessAnalysis
    %%% save out any important data
    postspikes{end+1:end+length(blockSpike)} = blockSpike;
    plot(blockSpike)


% for layer =2:6
% subplot(2,3,layer)
% bar(darkR)
% end