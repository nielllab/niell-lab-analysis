fname = 'allunits_new0108';

wv = xlsread(fname,1,'BU2:CM241');
inh = xlsread(fname,1,'P2:P241');
%use = find(~isnan(wv(:,1)));
%wv = wv(use,:);

inh_est = find((inh==1));
exc_est = find((inh~=1));

figure
hold on

for i=1:size(wv,1)

    if inh(i)==1;
        colorstr='r';
    elseif inh(i)==0;
        colorstr='g';
    else colorstr='.';
    end
    [y j] = min(wv(i,:));
    j
    if j <=5
        wvshift(i,:) = wv(i,1:17);
        plot(wv(i,1:17),colorstr);
    elseif j==6
           wvshift(i,:) = wv(i,2:18);
           plot(wv(i,2:18),colorstr);
    elseif j==7
        wvshift(i,:) = wv(i,3:19);
        plot(wv(i,3:19),colorstr);
    end
end


figure
imagesc(wvshift)

figure
hist(wvshift(:,11)./wvshift(:,17),[0:0.25:6])

figure
plot(inh(use),wvshift(:,11)./wvshift(:,17),'o')

[coeff scores latent] = princomp(wvshift);
figure
plot3(scores(:,1),scores(:,2),scores(:,3),'o')

[s coeff u] = fastica(wvshift(used,:)','numOfIC',2,'lastEig',2,'stabilization','on','g','tanh','approach','symm');
figure
plot(coeff);

figure
plot(s(1,find(inh(used)==1)),s(2,find(inh(used)==1)),'ro')
hold on
plot(s(1,find(inh(used)==0)),s(2,find(inh(used)==0)),'go')

figure
plot3(s(1,find(inh==1)),s(2,find(inh==1)),s(3,find(inh==1)),'ro')
hold on
plot3(s(1,find(inh==0)),s(2,find(inh==0)),s(3,find(inh==0)),'go')
% 
% figure
% plot(wv(inh_est,:)','r');
% hold on
% plot(wv(exc_est,:)','g');

