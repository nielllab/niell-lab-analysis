function xc = slidingXC(x,y,lagwin,slideWin)

pts = 1-slideWin(1):length(x)-slideWin(end);
clear xc
for i = 1:length(pts);
   i
   xc(i,:) = xcorr(x(pts(i)+slideWin),y(pts(i)+slideWin),'coeff',lagwin);
end
figure
imagesc(xc,[-0.5 0.5])
