function fit = fitProfile(sf,tf,ori,params,resp,shuffle)
if shuffle==1
    [shuffParams shuffResp] = resampleReps(params,resp,8);
else
    shuffParams=params; shuffResp=resp;
end

%%% perform shuffle
for s = 1:max(sf)
    for j=1:max(tf)
        for k=1:max(ori)
            trials = find(sf==s & tf==j & ori==k);
            respprofile(s,j,k,:)=shuffResp(trials);
        end
    end
end

%%% setup fit function
f = @(x)computeTuningErr(x,shuffParams,shuffResp);

%%% make initial guesses
meanprofile = squeeze(mean(respprofile,4));
sfguess = max(max(meanprofile,[],3),[],2) ;
tfguess = squeeze(max(max(meanprofile,[],3),[],1))';
origuess = squeeze(max(max(meanprofile,[],2),[],1));
x0 = [sfguess/max(sfguess); tfguess/max(tfguess); origuess; 0; 1] ;
%%% set upper lower bound
lb = zeros(18,1);lb = [lb ; -1; 0];
ub = [1 1 1 1 1 1 1 1 1 1 5 5 5 5 5 5 5 5 1 10 ]';

%%% do the fit ..
tic
fit = fmincon(f,x0,[],[],[],[],lb,ub);
toc
