figure
plot(movingsoff);
hold on
plot(movingcovered,'c');
plot(movingdim,'g');

plot(moving30percentbars,'m');
plot(movingbright,'k');
plot(movingsheet,'r');

figure
plot([movingsoff; movingcovered; movingdim; moving30percentbars; movingbright;movingsheet]')
legend('off','covered','dim','30%','full','sheet')