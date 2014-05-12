close all
use = find(OSI(:,1)>0.5 & OSI(:,2)>0.5 & peak(:,1)>3);

for i = 1:length(use)
    figure
    hold on
    for j = 1:2;
       t =  (squeeze(tuning(use(i),j,:)));
       t = [t(4:12); t(1:3)];
       plot(t)
    end
end

