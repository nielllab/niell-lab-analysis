%t = input('which subplot : ');
t=3;
subplot(3,3,t);

[x(1) x(2)] = (ginput(1));
x=round(x)
sta = get(gco,'CData');
x= max(x,16);
x= min(x,44);

subfield = sta(x(2)-15:x(2)+16,x(1)-15:x(1)+16);
figure
imagesc(subfield,[-64 64]);

