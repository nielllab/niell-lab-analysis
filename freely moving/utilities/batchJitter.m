%%% perform camera jitter analysis on a group of movie files
%%% loops through to select files, then runs checkCameraJitter on each

%%% select movie file names; click cancel to stop
done = 0;
n=0;
clear fname
while ~done
    n= n+1;
    [f p] = uigetfile('*.avi');
    if f~=0
        fname{n} = fullfile(p,f);
    else
        done = 1;
    end
end

%%% run checkCameraJitter
for i = 1:length(fname);
    i
    fname{i}
    [jitter(i,:) trace{i}] = checkCameraJitter(fname{i});
end

sprintf('x jitter = %0.2f +/- %0.2f pix',nanmean(jitter(:,1)), nanstd(jitter(:,1)))
sprintf('y jitter = %0.2f +/- %0.2f pix',nanmean(jitter(:,2)), nanstd(jitter(:,2)))

[f p] = uiputfile('*.mat','save results');
save(fullfile(p,f),'jitter','trace')
