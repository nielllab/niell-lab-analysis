%%%% run this section on each analysis file
selected =[];
for i = 1:size(contrast_tuning_centered,3);
    figure
    plot(contrast_tuning_centered(:,:,i));
    ginput(1);
    use = input('use this one? ');
    if use
       selected = [selected ; i];
    end
    contrast_tuning_centered(:,:,i) = contrast_tuning_centered(:,:,i)/max(max(contrast_tuning_centered(:,:,i)));
end
contrast_out = contrast_tuning_centered(:,:,selected);

save driftcon contrast_out
figure
plot(mean(contrast_out,3));


c_all = contrast_out;
c_all

%%%% then run this to select all the analysis files
contrastall_new=[];
done=0;

while done ==0
    [fname pname] = uigetfile('.mat');
    if fname ==0
        done=1;
    else
        load(fullfile(pname,fname),'contrast_out');
        sz = size(contrastall_new,3);
        sz2 = size(contrast_out,3);
        if isempty(contrastall_new)
            contrastall_new=contrast_out;
        else
            if size(contrast_out,2)==5
                contrastall_new(:,1:5,sz+1:sz+sz2) = contrast_out;
            else
                 contrastall_new(:,2:5,sz+1:sz+sz2) = contrast_out;
                 contrastall_new(:,1,sz+1:sz+sz2) =NaN;
            end
        end
    end
end

% %%% 060707 unit 5_1 is the exception
% contrastall(:,:,14)=NaN;

%%% run this to generate figures

max(contrastall(:,:,9))



for i = 1:size(contrastall,3);
   for i = 17
    for j = 1:5

        theta_ind = [0:11]*pi/6;
        theta_tuning = squeeze(contrastall(:,j,i))';
        if isnan(theta_tuning)
            peak(i,j)=NaN;
            theta(i,j) = NaN;
            osi(i,j)=NaN;
            width(i,j)=NaN;
            A1(i,j)=NaN;
            B(i,j)=NaN;
        else
            peak(i,j) = max(theta_tuning);
            [theta(i,j) osi(i,j) A1(i,j) A2 ...
                width(i,j) B(i,j) null(i,j) yfit]= fit_tuningcurve(theta_tuning,theta_ind);
        end
    end
end
null_original=null;
width_original = width;

width = width*180/pi;
width = abs(width);

peak = A1+B;
null(null<-.05)=0;
osi = (peak-null)./(peak+null);
use = peak(:,5)>0.5 & peak(:,5)<2;
use(8)=0;


figure
plot(width(use,4),width(use,5),'o');
hold on
plot([0 40],[0 40]);
figure
plot(osi(use,4),osi(use,5),'o');
hold on
plot([0 1.2],[0 1.2]);
mean(osi(use,:))
std(osi(use,:))/sqrt(32)

figure
plot(A1'+B');
axis([1 5 0 1.1])





peak = A1+B;

for i = 1:35;
   peak(i,:)= peak(i,:)./max(peak(i,:));
end

figure
plot(peak');

figure
plot(peak(use,:)')
cr = zeros(6,sum(use));
cr(2:6,:) = peak(use,:)';
figure
plot([0 6 12 25 50 100],cr);

figure
plot(peak')

used = (osi(:,5)>0.5 & osi(:,4)>0.5) &use;
width(width>45)=45;
figure
plot(width(used,4),width(used,5),'o');
hold on
plot([0 50],[0 50])

mean(width(used,:))
std(width(used,:))/sqrt(sum(used))

contrastplot = contrastall(:,2:5,used);
figure
plot(mean(contrastplot,3))
for i = find(used==1)';
   figure
   plot(contrastall(:,:,i))
   title(sprintf('used %d',i));
end


newwidth = sqrt(-2*log(sqrt(0.5) - (1-osi)./(1+osi))).*width;
newwidth(~isreal(newwidth))=NaN;
newwidth(newwidth>60)=60;
figure
plot(newwidth(used,4),newwidth(used,5),'o');
mean(abs(newwidth(used,:)))
std(abs(newwidth(used,:)))/sqrt(sum(used))

hold on
plot([0 70],[0 70])

for i = 1:35;
   peak(i,:)= peak(i,:)./max(peak(i,:));
end
use = peak(:,5)>0.5 & peak(:,5)<2;
use(8)=0;
cr = peak(:,:);
contrasts = [6.25 12.5 25 50 100];
for i = 1:35;
 
    ctuning = cr(i,:);
    upper = find(ctuning>0.5);
    p0(1)=contrasts(min(upper));
    if p0(1)>50;
        p0(1)=50;
    end
    p0(2)= 1;
    p = nlinfit(contrasts,ctuning,@naka_rushton,p0);
    halfC(i) = p(1);
    slope(i) = p(2);
    figure
    plot(ctuning);
    hold on
    plot(naka_rushton(p,contrasts),'g');
end
    
figure
hist(slope(use));

figure
hist(log2(halfC(use)/6.25),0:0.5:3);



halfslope = 100*0.25*slope(use)./halfC(use);
halfslope(halfslope>8)=8;
figure
hist(halfslope,0.5:1:10);


figure
plot(slope(use),halfC(use),'o')

figure
plot(contrasts,cr')

figure
lowslope = find(slope<2);
plot(contrasts,cr(lowslope,:)');

figure
for i = 1:size(cr,1);
%     figure
%     plot(contrasts,cr(i,:)');
   if use(i)
       hold on
     plot(contrasts,naka_rushton([halfC(i) slope(i)],contrasts),'g')
    title(sprintf('%d',i));
   end
end

plotcr = [  31 27]
figure
for i = plotcr
    hold on
    plot(contrasts,cr(i,:)','o');
    plot(contrasts,naka_rushton([halfC(i) slope(i)],contrasts),'g')
end    

figure
cr_all= zeros(6,sum(use))
cr_all(2:6,:) = cr';
plot(cr_all)



s=~isnan(contrastall);
contrastall(isnan(contrastall))=0;
cmean = sum(contrastall,3)./sum(s,3);
cstd = sqrt(sum(contrastall.^2,3) - sum(s,3).*cmean.^2)./sum(s,3)
figure
plot(cmean)
figure
errorbar(cmean,cstd)

cmean = mean(contrastplot,3);
cstd = std(contrastplot,[],3)./sqrt(sum(used));
cmeanshift(1:3,:)=cmean(10:12,:);
cmeanshift(4:12,:) = cmean(1:9,:);
cstdshift(1:3,:)=cstd(10:12,:);
cstdshift(4:12,:) = cstd(1:9,:);
figure
plot(cmeanshift(:,1:4));
legend('12%','25%','50%','100%')

figure
errorbar(cmeanshift(:,2:5),cstdshift(:,2:5));
legend('12%','25%','50%','100%')

for i=1:5;
    m = max(cmean(:,i)-min(cmean(:,i)));
    cmeannorm(:,i) = (cmeanshift(:,i)-min(cmean(:,i)))/m;
    cstdnorm(:,i) = cstd(:,i)/m;
end
figure
plot(cmeannorm(:,2:5))
legend('12%','25%','50%','100%')

figure
errorbar(cmeannorm(:,2:5),cstdnorm(:,2:5));
legend('12%','25%','50%','100%')

for i = 1:size(contrastall,3);
    figure
    plot(contrastall(:,:,i));
    title(sprintf('unit %d',i));
end
    