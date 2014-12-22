close all
clear all

%%% make fake data  - list of R, sf, tf, ori, sp
dF=0.05;
sftuning = [1 2 5 10 2]; sftuning = sftuning/max(sftuning);
tftuning = [1 2 6 5 1]; tftuning = tftuning/max(tftuning);
osituning = dF*[1 3 10 2 1 5 20 2]/4;
gain=1;
nreps = 8;

tr=0;
for i = 1:length(sftuning)
    for j=1:length(tftuning)
        for k=1:length(osituning)
            for r = 1:nreps
                tr = tr+1;
                sp(tr) = round(rand);
                R(tr) = sftuning(i)*tftuning(j)*osituning(k)*(1+gain*sp(tr));
                sf(tr)=i; tf(tr)=j; ori(tr)=k;
            end
        end
    end
end

%%% model as poisson with noise
R = dF*poissrnd(R/dF);
R = R+normrnd(0,2*dF,size(R));


%%% set up data - move parameters (ori/sf/tf/run) into one array, 
%%% put responses in matching order

tr=1;
for i = 1:length(sftuning)
    for j=1:length(tftuning)
        for k=1:length(osituning) 
                params(tr:tr+nreps-1,i)=1;params(tr:tr+nreps-1,j+5)=1; params(tr:tr+nreps-1,k+10)=1;            
                trials = find(sf==i & tf==j & ori==k);
                params(tr:tr+nreps-1,19)=sp(trials);
                resp(tr:tr+nreps-1) = R(trials);
            respprofile(i,j,k,:)=R(trials);
            tr = tr+nreps;
        end
    end
end


figure
imagesc(params)

%%% show sf/tf tuning for 8 orientations
layout = [6 3 2 1 4 7 8 9];
meanprofile = squeeze(mean(respprofile,4));
figure
range = prctile(meanprofile(:),99)
if range<dF
    range=dF;
end
for o = 1:8
    subplot(3,3,layout(o))
    imagesc(squeeze(meanprofile(:,:,o)),[0 range])
    axis off
    axis square
    axis equal
end


gcp  %%% opens matlab pool
clear sfguess origuess tfguess
tic
N=100
parfor i = 1:N %%% repeat N times for bootstrap
i
if i==1 %%% first time is not shuffeld
    shuffle=0;
else shuffle=1;
end

fit = fitProfile(sf,tf,ori,params,resp, shuffle); %%% performs the fit!

%%% scale amplitude of sf/tf so responsiveness is all in ori
amp1 = max(fit(1:5)); amp2 = max(fit(6:10));
sfest(i,:) = fit(1:5)/amp1;
tfest(i,:) = fit(6:10)/amp2;
oriest(i,:) = fit(11:18)*amp1*amp2;

%%% get peak sf and tf, and max of ori (which is responsiveness)
[y ind] = max(sfest(i,:));
sfmax(i) = ind;
[y ind] = max(tfest(i,:));
tfmax(i)=ind;
orimax(i) = max(oriest(i,:));
[allOSI(i) alltheta(i)]= calcOSI(oriest(i,:)',0);
base(i) = fit(19);  
allGain(i) = fit(20);
end
toc

sprintf(' sf = %f +/- %f',mean(sfmax),std(sfmax))
sprintf(' tf = %f +/- %f',mean(tfmax),std(tfmax))
sprintf(' ori = %f +/- %f',mean(allOSI),std(allOSI))
sprintf('resp = %f +/- %f',mean(orimax+base),std(orimax))
sprintf('p < %fof resp less than 0.1',sum((orimax+base)<0.1)/length(orimax))
sprintf('theta = %f +/- %f',mean(alltheta)*180/pi,std(alltheta)*180/pi)
sprintf('base = %f +/- %f',mean(base),std(base))
sprintf('running gain = %f +/- %f; actual %f',mean(allGain),std(allGain), gain)

figure
errorbar(sfest(1,:),std(sfest,1))
hold on
plot(sftuning,'r')
legend('fit','actual')
xlabel('SF')

figure
errorbar(tfest(1,:),std(tfest,1))
hold on
plot(tftuning,'r')
legend('fit','actual')
xlabel('TF')

figure
errorbar(oriest(1,:)+base(1),std(oriest,1))
xlabel('orientation')
hold on
plot(osituning,'r')
legend('fit','actual')

% figure
% hist(allOSI)
% xlabel('OSI'); xlim([0 1])
% ylabel('N')


