function PlotSpontData(data1,data2,treatment,lyr,inh,resp,titlestr)

figure
plot(data1(treatment &lyr & resp),data2(treatment & lyr & resp),'o')
hold on
plot(data1(treatment & lyr &  resp & inh),data2(treatment& lyr & resp & inh),'ro','Linewidth',2)
plot(data1(treatment & lyr & resp &~inh),data2(treatment & lyr & resp & ~inh),'go','Linewidth',2)
hold on
plot([0 40],[0 40])
axis equal
legend('all','inh','exc')
title (titlestr)
xlabel('pre'); ylabel('post')