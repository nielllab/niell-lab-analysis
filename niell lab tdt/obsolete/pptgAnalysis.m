% pname = uigetdir('C:\data\TDT tanks','block data')
% delims = strfind(pname,'\');
% selected_path = pname(1 :delims(length(delims))-1)
% Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
% Block_Name = pname(delims(length(delims))+1 :length(pname))

 Tank_Name='06222011_mlr_stim_awake'
 Block_Name='wn3_72hz'
chans =1:4:16;
flags = struct('lfpTseries',1,'lfpSpectra',1','mouseOn',1,'laserOn',1,'MUspike',1,'visStim',1)

tdtData= getTDTdata(Tank_Name, Block_Name, chans, flags);
%tdtData= getTDTdata(Tank_Name, Block_Name, nChan, 0, 0, 0, 1,1, 0, 0, 0);


smoothwindow_secs = 0.1;
laserT = tdtData.laserT;
dt = laserT(2)-laserT(1);
smoothwindow = ceil(smoothwindow_secs/dt)
lasersmooth = zeros(size(tdtData.laserTTL));
%%% replace this with convolution (for speed)!
for i = smoothwindow+1:length(tdtData.laserTTL);
    lasersmooth(i)= mean(tdtData.laserTTL(i-smoothwindow:i));
end
lasersmooth = lasersmooth/5;
figure
plot(lasersmooth)

tic
tsamp = 0;
npulse = 0;
for i = 2:length(lasersmooth);
    if lasersmooth(i-1)==0 & lasersmooth(i)>0
        npulse = npulse+1;
        onT(npulse) = laserT(i);
    end
end
toc

timeRange = -10:0.1:25;
vdata=zeros(length(timeRange),npulse-1);
for i = 2:npulse-2
    vdata(:,i) = interp1(tdtData.mouseT,tdtData.mouseV,timeRange+onT(i));
end
figure
plot(timeRange,vdata)
hold on
plot([0 17],[40 40],'g','LineWidth',12)

figure
plot(timeRange,mean(vdata,2))
hold on
plot([0 17],[20 20],'g','LineWidth',12)


for ch = chans;
    df = median(diff(tdtData.spectF{ch}));
    specdata = tdtData.spectData{ch};
    tdata = tdtData.spectT{ch};
    spect_avg = zeros(length(timeRange),size(specdata,2));
    for i = 2:npulse-1;
        for f=1:size(specdata,2);
            spect_avg(:,f) = spect_avg(:,f) + interp1(tdata,squeeze(specdata(:,f)),timeRange+onT(i))';
        end
    end
    figure
    imagesc(spect_avg',[0 2*10^-9])
    axis xy
     df = median(diff(tdtData.spectF{ch}));
     dt = median(diff(timeRange));
     set(gca,'YTick',(10:10:80)/df);
        set(gca,'YTickLabel',{'10','20','30','40','50','60','70','80'})
         set(gca,'XTick',(5:5:40)/dt);
        set(gca,'XTickLabel',{'-5','0','5','10','15','20','25','30'})
end


bins = 0:0.25:max(tdtData.mouseT);
R= zeros(length(bins),nChan);
for ch = chans;
    R(:,ch)= hist(tdtData.MUspikeT{ch},bins);
end
figure
plot(R)

f= tdtData.frameEpocs(1,:);
t= tdtData.frameEpocs(2,:);
laserInterp = interp1(laserT,lasersmooth,t);
figure
hist(diff(t))
onT
dt = (diff(t));
dt_avg = median(dt);
for ch = chans
    for rep = 1:2
        R=zeros(1,max(f));
        cycframes = 300;
        Rcycle = zeros(1,cycframes);
        ntrial = zeros(1,max(f));
        cyctrials = zeros(1,cycframes);
        
        for i = 1:length(tdtData.frameEpocs)-1
           i
           if rep ==1 & laserInterp(i)==0
                ntrial(f(i))=ntrial(f(i))+1;
                R(f(i)) = R(f(i))+sum(tdtData.MUspikeT{ch}>t(i) & tdtData.MUspikeT{ch}<t(i)+dt_avg);
                Rcycle(mod(f(i),cycframes)+1) = Rcycle(mod(f(i),cycframes)+1) + sum(tdtData.MUspikeT{ch}>t(i) & tdtData.MUspikeT{ch}<t(i)+dt_avg);
                cyctrials(mod(f(i),cycframes)+1) = cyctrials(mod(f(i),cycframes)+1) +1;
            elseif rep==2 & laserInterp(i)>0.3;
                ntrial(f(i))=ntrial(f(i))+1;
                R(f(i)) = R(f(i))+sum(tdtData.MUspikeT{ch}>t(i) & tdtData.MUspikeT{ch}<t(i)+dt_avg);
                Rcycle(mod(f(i),cycframes)+1) = Rcycle(mod(f(i),cycframes)+1) + sum(tdtData.MUspikeT{ch}>t(i) & tdtData.MUspikeT{ch}<t(i)+dt_avg);
                cyctrials(mod(f(i),cycframes)+1) = cyctrials(mod(f(i),cycframes)+1) +1;
            end
        end
        R = R./ntrial;
        figure
        plot(R)
        Rcyclerep{rep} = Rcycle./cyctrials;
        figure
        plot(Rcycle)
        %    figure
        %    plot(1-cos(2*pi*(1:300)/300),Rcycle)
    end
    figure
    plot(Rcyclerep{1});
    hold on
    plot(Rcyclerep{2},'g');
end






