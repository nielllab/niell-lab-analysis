close all
clear correct
clear correctOrient
data =[];
for fitsf=1:2
    nreps=8;
    
    for ori = 1:8
        data = [data ; squeeze(allresp(:,ori,fitsf,:))'];
        %orientation((ori-1)*nreps + (1:nreps))=mod((ori-1),4)+1;
    end
end
sf(1:64)=1; sf(65:128) =2;

figure
imagesc(data)
figure
plot(sf)
for sz =1:5
    sampSize = 4^(sz-1);
    if sampSize>size(data,2)
        sampSize = size(data,2)-1
    end
    sampSize
    tic
    counts = zeros(8,1); testcount = zeros(size(data,1),1); trialcorrect=testcount;
    confusion=zeros(8,8);
    for szIter = 1:40
        useSamps = randsample(size(data,2),sampSize);
        useData = data(:,useSamps);
        
        nfold=10;
        c = cvpartition(size(data,1),'k',nfold);
        
        for iter = 1:nfold
            
            sv = fitctree(useData(c.training(iter),:),sf(c.training(iter)));
            %sv = fitensemble(useData(c.training(iter),:),orientation(c.training(iter)),'Subspace',64,'discriminant');
            label = predict(sv,useData(c.test(iter),:));
            shouldlabel=sf(c.test(iter));
            tr=find(c.test(iter));
            for l = 1:length(label);
                counts(shouldlabel(l))=counts(shouldlabel(l))+1;
                confusion(shouldlabel(l),label(l))= confusion(shouldlabel(l),label(l))+1;
                testcount(tr(l))=testcount(tr(l))+1;
                if label(l)==shouldlabel(l)
                    trialcorrect(tr(l))=trialcorrect(tr(l))+1;
                end
            end
            
            correct(sz,szIter,fitsf,iter) = sum(label' == sf(c.test(iter)))/length(label);
            correctOrient(sz,szIter,fitsf,iter) = sum(mod(label',4) == mod(sf(c.test(iter)),4))/length(label);
        end
        
    end
    toc
    figure
    repcounts = repmat(counts',8,1);
    imagesc(confusion./repcounts,[0 1]);
    figure
    plot(trialcorrect./testcount)
    ylim([0 1])
end


avgCorrect = squeeze(mean(mean(correct,4),2))
avgCorrectOrient = squeeze(mean(mean(correctOrient,4),2))

stdCorrect = squeeze(std(mean(correct,4),[],2))
stdCorrectOrient = squeeze(std(mean(correctOrient,4),[],2))
figure
errorbar([1 2 3 4 4.5],avgCorrect(:,1),stdCorrect(:,1)/sqrt(size(correct,2)))
hold on
errorbar([1 2 3 4 4.5],avgCorrect(:,2),stdCorrect(:,2)/sqrt(size(correct,2)),'g')
ylim([0 1])

figure
errorbar([1 2 3 4 4.5],avgCorrectOrient(:,1),stdCorrectOrient(:,1)/sqrt(size(correct,2)))
hold on
errorbar([1 2 3 4 4.5],avgCorrectOrient(:,2),stdCorrectOrient(:,2)/sqrt(size(correct,2)),'g')
ylim([0 1])

