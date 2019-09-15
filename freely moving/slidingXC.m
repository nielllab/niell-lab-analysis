function xc = slidingXC(x,y,lagwin,slideWin)

pts = 1-slideWin(1):length(x)-slideWin(end);
clear xc
for i = 1:length(pts);
   xc(i,:) = nanxcorr(x(pts(i)+slideWin),y(pts(i)+slideWin),lagwin,'zero');
end

