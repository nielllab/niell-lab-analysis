%jumpy_plotting.m
load('D:\labeledDLC\pdfs\CK2-CK2-7P1-RT_AllSessions_112719.mat')

for i = 1:size(Data,2)
    figure
    plot(Data(i).Rtheta-nanmean(Data(i).Rtheta))
    hold on
    plot(Data(i).Ltheta-nanmean(Data(i).Ltheta))
    plot(Data(i).theta*(180/pi)-nanmean(Data(i).theta*(180/pi)))
    legend('Reye theta','Leye theta','head theta')
end

close all
for i = 1:size(Data,2)
    figure
    plot(Data(i).XRcent-nanmean(Data(i).XRcent))
    hold on
    plot(Data(i).XLcent-nanmean(Data(i).XLcent))
    plot(Data(i).theta*(180/pi)-nanmean(Data(i).theta*(180/pi)))
    legend('Reye theta','Leye theta','head theta')
end