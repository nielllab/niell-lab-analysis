%%% test approach
%load('Analyzed_AllAnimals_090319_New.mat')
trialNum=20  %this trial is not stitching correctly

a=ani{useData(trialNum)}; d=date{useData(trialNum)};s=sess{useData(trialNum)};c=clipnum(useData(trialNum));
disp([a,d,s,c])

appT=find(appEpoch{trialNum}==1); startT=appT(1)/30

figure;plot(mouseSp{useData(trialNum)}); hold on;
plot(appT,mouseSp{useData(trialNum)}(appEpoch{trialNum}==1),'g')