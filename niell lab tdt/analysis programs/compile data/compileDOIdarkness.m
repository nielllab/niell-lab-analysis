clear all
close all

dbstop if error

batchDOIephys; %%% load batch file

%%% select the sessions you want based on filters
%%% example
use =  find(strcmp({files.treatment},'Saline') & strcmp({files.notes},'good data') & ~cellfun(@isempty,{files.predark})   )

%use =  find( strcmp({files.treatment},'DOI') & strcmp({files.notes},'good data') & ~cellfun(@isempty,{files.predark}) & ~cellfun(@isempty,{files.postdark}) )

sprintf('%d selected sessions',length(use))

cv2all=[]; %%% initialize data to be saved
meanRwn = []; meanRdark = [];
preCorrWN = []; postCorrWN=[]; preCorrDark = []; postCorrDark=[];
layers = [];

savePDF=0;
redo = 1;
for i = 1:2
    afile = [pathname '\' files(use(i)).dir '\' files(use(i)).analysisfile '.mat'];
    clustfile = [pathname '\' files(use(i)).dir '\' files(use(i)).clusterfile '.mat'] ;
     
    [inh mid] = getWaveform(clustfile,afile,0);
    
    lfpMove = getLfpMovement(clustfile,afile,files(use(i)).blockWn{1},0);
    preLFP(i,:,:) =squeeze(median(lfpMove.meanSpect, 1))/median(lfpMove.meanSpect(:));
    
    lfpMove = getLfpMovement(clustfile,afile,files(use(i)).blockWn{2},0);
    postLFP(i,:,:) =squeeze(median(lfpMove.meanSpect, 1))/median(lfpMove.meanSpect(:));
    
    dt = 1;
    [preCorr postCorr cv2 R eigs] = prepostDOIdarkness(clustfile,afile,files(use(i)).blockWn,dt,0);
    preCorrWN = [preCorrWN; preCorr(:)]; postCorrWN = [postCorrWN; postCorr(:)];
    meanRwn = [meanRwn ; squeeze(mean(R,2))]
    
    
    dt = 1;
    [preCorr postCorr cv2 R eigs] = prepostDOIdarkness(clustfile,afile,files(use(i)).blockDark,dt,0);
    preCorrDark = [preCorrDark; preCorr(:)]; postCorrDark = [postCorrDark; postCorr(:)];
    
    cv2all = [cv2all; cv2];
    meanRdark = [meanRdark ; squeeze(mean(R,2))];
    
    eigsFull = zeros(100,2);
    eigsFull(1:size(eigs,1),:)=eigs;
    eigsAll(i,:,:) = eigsFull;
   
    load(afile,'layer');
    layers = [layers; layer];
end

figure; hold on
for i = 1:size(preLFP,1);
    plot(preLFP(i,:,1),'r');
    plot(preLFP(i,:,2),'g');
end
title('pre')

figure; hold on
for i = 1:size(preLFP,1);
    plot(postLFP(i,:,1),'r');
    plot(postLFP(i,:,2),'g');
end
title('post')

