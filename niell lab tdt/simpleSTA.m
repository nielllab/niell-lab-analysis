%%% calculate STA from spike data, with 9 temporal lags
%%% needs:
%%% spikehist(frame) = # of spikes fired in each from of the movie
%%% moviedata(x,y,frame) = movie image frames

%%%
sta_all= zeros(size(moviedata,1),size(moviedata,2),9);  %%% spatiotemporal STA

lagframes= find(spikehist>0); %%% to save time, only perform calculations on frames with at least one spike
lagframes=lagframes(lagframes>=10);   %%% can't do STA on first 9 frames since don't have enought lags
for lag = 1:9  %%% loop over each timelag
    for i = lagframes  %%% loop over all the frames with spikes
        %%% calculate average by weighted sum
        sta_all(:,:,lag) =sta_all(:,:,lag)+spikehist(i)*moviedata(:,:,fnum(i-lag));
    end
end

%%% normalize and zero-mean
sta_all = sta_all/sum(spikehist);
sta_all = (sta_all-128)/128;

%%% perform SVD to separate into spatial and temporal components
%%% (a simpler option if you only want spatial component is to just choose one timelag)

%%% reshape into space,time rather than x,y,time
obs = reshape(sta_all,size(sta_all,1)*size(sta_all,2),size(sta_all,3));
[u s v] = svd(obs');

%%% select the first spatial and temporal component
meanSTA= reshape(v(:,1),size(sta_all,1),size(sta_all,2));
tcourse = u(:,1);
if sum(tcourse(1:5))<0;
    tcourse = -tcourse;
    meanSTA=-meanSTA;
end
meanSTA = meanSTA*s(1,1)*max(tcourse);
tcourse = tcourse/max(tcourse);

%%% display the results
subplot(1,2,1);
imagesc(meanSTA,[-0.1 0.1]); axis square; axis off; colormap jet
subplot(1,2,2);
plot(tcourse); hold on; plot([1 9],[0 0],':'); ylim([-1.05 1.05])