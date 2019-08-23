load('J463bAllVids.mat')
win = -5:5;
clear x y
i=1
 
x = Data(i).dxLTheta;
y = Data(i).dxRTheta;
win = -5:5;
bad = abs(Data(i).Rtheta)>45 | abs(Data(i).Rphi)>45; bad = bad(1:length(x));
[b yfit r] = fitEyeData(x,y,win,bad);
title('R dTh vs L dTh')

x = Data(i).dxLPhi;
y = Data(i).dxLTheta;
win = -5:5;
bad = abs(Data(i).Ltheta)>45 | abs(Data(i).Lphi)>45; bad = bad(1:length(x));
[b yfit r] = fitEyeData(x,y,win,bad);
title('dPhi vs dTh')

x = Data(i).dtheta(1:end-1)*180/pi;
yfit = [0; 0; 0; 0; 0; yfit; 0; 0; 0; 0; 0];
y = y - yfit;
win = -5:5;
bad = abs(Data(i).Rtheta)>45 | abs(Data(i).Rphi)>45; bad = bad(1:length(x));
[b yfit r] = fitEyeData(x,y,win,bad);
title('fit to dth after dphi is removed')

x = Data(i).dtheta(1:end-1)*180/pi;
y = Data(i).dxLPhi;
win = -5:5;
bad = abs(Data(i).Rtheta)>45 | abs(Data(i).Rphi)>45; bad = bad(1:length(x));
[b yfit r] = fitEyeData(x,y,win,bad);
title('head dth vs R dphi')

x = Data(i).dtheta(1:end-1)*180/pi;
y = Data(i).dxRTheta;
win = -5:5;
bad = abs(Data(i).Rtheta)>45 | abs(Data(i).Rphi)>45; bad = bad(1:length(x));
[b yfit r] = fitEyeData(x,y,win,bad);
title('head dth vs R dth')

angleRange = -90:10:90;
for j = 1:length(angleRange);
    x = Data(i).dtheta(1:end-1)*180/pi;
    y = Data(i).dxRTheta*cosd(angleRange(j)) + Data(i).dxRPhi*sind(angleRange(j));
    win = -5:5;
    bad = abs(Data(i).Rtheta)>45 | abs(Data(i).Rphi)>45; bad = bad(1:length(x));
    [b yfit r(j)] = fitEyeData(x,y,win,bad);
    title(sprintf('rotation angle = %d',angleRange(j)))
end
figure
plot(angleRange,r)


 bad = abs(Data(i).Rtheta)>45 | abs(Data(i).Rphi)>45 | abs(Data(i).Ltheta)>45 | abs(Data(i).Lphi)>45; bad = bad(1:end-1);

x= [Data(i).dxRTheta(~bad)  Data(i).dxLTheta(~bad)   Data(i).dxRPhi(~bad)  Data(i).dtheta(~bad)*180/pi;];
y= Data(i).dxLPhi(~bad);

b = glmfit(x,y,'normal')
yfit = glmval(b,x,'identity');
r = sqrt(1-nansum((y-yfit).^2)/nansum(y.^2))

lagwin = 10;
slideWin = -30:30;

x = Data(i).dtheta(1:end-1)*180/pi;
y = Data(i).dxRTheta;
xc = slidingXC(x,y,lagwin,slideWin);

x = Data(i).dtheta(1:end-1)*180/pi;
y = Data(i).dxRPhi;
xc = slidingXC(x,y,lagwin,slideWin);

x = Data(i).dxRTheta;
y = Data(i).dxRPhi;
xc = slidingXC(x,y,lagwin,slideWin);

x = Data(i).dxRTheta;
y = -Data(i).dxLTheta;
xc = slidingXC(x,y,lagwin,slideWin);


