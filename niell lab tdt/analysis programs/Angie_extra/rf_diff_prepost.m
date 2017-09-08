
%%
% figure
% for i = 1:length(data_wn)
% x = abs(diff(sta_all_img{i,1}-sta_all_img{i,2}));
% subplot(3,10,i)
% imagesc(x); axis square
% % x = reshape(x,size(x,1)*size(x,2),size(x,3));
% % x =(mean(x,1));
% end


%%


clear max_wn data_wn
% low_wn = squeeze(min(wn_crf,[],2));
%amp_wn = wn_evoked-wn_spont; 
amp_wn = evoked-spont;
ampmean = squeeze(mean(amp_wn,2))
useResp_wn = amp_wn(:,2,1)>3| amp_wn(:,2,2)>3; %| amp_wn(:,2,1)>3 | amp_wn(:,2,2)>3;
%useResp_wn = ampmean(:,1)>3 | ampmean(:,2)>3
data_wn = goodAll==1 &hasWn==1 & useResp_wn' & inhAll==0;


% cc = zeros(length(find(data_wn==1 & treatment==1 & hasWn==1)),121,121);

titles={'Saline','DOI'};
clear cc_all
figure
for t= 1:2
use = find(data_wn==1 &treatment==t);

% cc = zeros(length(use),121,121);
 %ccmean= zeros(length(use));
for i = 1:length(use)
a=sta_all_img{use(i),1}; a = a(:);
b=sta_all_img{use(i),2}; b = b(:);
%cc(i,:,:) = xcorr2(a,b);
 [h1(i) p(i)]= ttest(a,b);
 cc = corrcoef(a,b); cc = cc(2,1);
 cc_all(t,i) = cc;
% ccmean= mean(mean(cc_all,2),3);
%af = fft(a); bf=fft(b)
% subplot(3,5,i);plot(af); %ylim([-50 50]); 
% hold on;
% plot(bf,'r'); %ylim([-50 50]);
%plot(a);% hold on;plot(b) %ylim([-50 50]); 

end
subplot(3,2,t)
ccmean = squeeze(mean(mean(cc,2),3));
Mbins=-1:.1:1;
h = hist(cc_all(t,:),Mbins);axis square; %ylim([0 15]);
miWn= (mean(wn_evoked(use,2,2),3)-mean(wn_evoked(use,2,1),3))./...
    (mean(wn_evoked(use,2,2),3)+mean(wn_evoked(use,2,1),3));
% h= hist(miWn,-1:.2:1); 
title(titles{t});
bar(Mbins,h/length(use),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .3]);
ylabel('proportion of cells');

subplot(3,2,t+2)
nh = hist(h1,0:1:1); bar(0:1:1,nh/length(h1)); axis square;% ylim([0 1]);
ylabel('proportion of cells'); %0 = reject null if no difference

 subplot(3,2,t+4);
 %xlim([-1.5 1.5]);
% 
scatter(miWn,cc_all(t,:)); ylabel('cc value');xlabel('MI');axis square;
xlim([-1.5 1.5]);
% xlim([-6 10]);
end

%%

%figure; subplot(1,2,1);plot(a); ylim([-50 50]);subplot(1,2,2); plot(b); ylim([-50 50]);
% usea =find( a>(.2*(range(a))))% | a<(.2*(range(a)))
% figure;plot(a(usea))
% useb = b<(.2*(range(b)))% | a<(.2*(range(a)))
% figure;plot(b(useb))
figure
for t = 1:2
use = find(data_wn==1 &treatment==t);

for i = 1:length(use)
a=sta_all_img{use(i),1}; a = a(:);
b=sta_all_img{use(i),2}; b = b(:);
qa = quantile(a(i),[0.025 0.25 0.50 0.75 0.975]);
qb = quantile(b(i),[0.025 0.25 0.50 0.75 0.975]);


usea =find( a<(qa(2)) | a>(qa(4)));
useb =find( a<(qb(2)) | a>(qb(4)));


% figure;
% subplot(121); plot(a(usea)); axis square; ylim([-50 50]); 
% subplot(122); plot(b(useb));axis square; ylim([-50 50]); 
% 
% 
end
 c = corrcoef(a(usea),b(useb(1:length(usea)))); c = c(2,1);
 cc_all(t,i) = cc;


subplot(1,2,t)
 Mbins=-1:.1:1;
 h = hist(cc_all(t,:),Mbins);axis square; %ylim([0 15]);
bar(Mbins,h/length(use),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .5]); axis square

end
%test = anova1(cc_all(1,:),cc_all(2,:))
%%
%MI for wn evoked
for t=1:2
  
    useN = data_wn & treatment==t
    miWn= (mean(wn_evoked(useN,2,2),3)-mean(wn_evoked(useN,2,1),3))./...
        (mean(wn_evoked(useN,2,2),3)+mean(wn_evoked(useN,2,1),3));
    h= hist(miWn,-1:.2:1);
    Mbins=-1:.2:1;
    subplot(1,2,t);
    bar(Mbins,h/sum(useN),'FaceColor',[0 .5 .5],'Linewidth',2);%ylim([0 .4]);
    xlim([-1.5 1.5]);axis square; 
    end


%%
% for change in FR, what is cc?
figure
scatter(miWn,ccmean)
%ylim([-.5 .5]);
%%
figure
clear r
for i = 1:length(data_wn)
subplot(3,10,i);
imagesc(sta_all_img{i,1});colormap jet
axis square
end

figure
clear r
for i = 1:length(data_wn)
subplot(3,10,i);
imagesc(sta_all_img{i,2});colormap jet
axis square
end

%%%%fit STAs%%%%
%%

clear r
figure
for i = 1:length(data_wn)
r = xcorr2(sta_all_fit{i,1},sta_all_fit{i,2});
subplot(3,10,i)
imagesc(r); axis square
end

figure
clear r
for i = 1:length(data_wn)
subplot(3,10,i);
imagesc(sta_all_fit{i,1})
axis square
end

figure
clear r
for i = 1:length(data_wn)
subplot(3,10,i);
imagesc(sta_all_fit{i,2})
axis square
end
