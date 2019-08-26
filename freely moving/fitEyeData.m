function [b yfit r] = fitEyeData(xIn,yIn,win,bad)

y=yIn; y(bad)=NaN;

pts = 1-win(1): length(xIn)-win(end);
y = y(pts);

for j= 1:length(pts)
    x(j,:) = xIn(pts(j)+win);
end

b = glmfit(x,y,'normal')
yfit = glmval(b,x,'identity');
figure
subplot(2,2,1)
plot(xIn(~bad),yIn(~bad),'.'); xlabel('x'),ylabel('y'); axis equal
hold on
%plot(xIn(~bad),glmval(b,xIn(~bad),'identity'),'.');

subplot(2,2,2);
plot(win,b(2:end)); hold on; plot(win,zeros(size(win)),'r:');
r = sqrt(1 - nansum((y-yfit).^2)/nansum(y.^2))
title(sprintf('b-weights r = %0.3f',r));

subplot(2,2,3)
[xc lag] = xcorr(xIn(~bad),yIn(~bad),'coeff',max(abs(win)));
plot(lag,xc); hold on; plot(lag,zeros(size(lag)),'r:');
title('xcorr');

subplot(2,2,4)
    plot(y,yfit,'.'); axis equal; hold on;plot([-10 10],[-10 10],'Linewidth',2)
    xlabel('y'); ylabel('yfit')

