function plotPinpData(data,wt,lyr,pinped,inh,resp)

figure
plot(data(wt &lyr & resp,1),data(wt & lyr& resp,2),'o')
hold on
plot(data(wt & lyr &  resp & inh,1),data(wt& lyr & resp & inh,2),'mo')
plot(data(wt & lyr & pinped & resp &~inh,1),data(wt & pinped &lyr & resp& ~inh,2),'go','Linewidth',2)
plot(data(wt & lyr & pinped & resp & inh,1),data(wt & pinped & lyr & resp & inh,2),'ro','Linewidth',2)
axis equal
legend('exc','inh','pinp exc','pinp inh')
xlabel('stop'); ylabel('run')