function z = gauss2d(p,x);
%%% calculate 2d gaussian based on parameters
%%% p = [A B x0 y0 sigx sigy], where ...
%%% A = amplitude of baseline
%%% B = baseline
%%% x0,y0 = center coordinates
%%% sigx, sigy = width of gaussian

a= p(1);
b=p(2);
x0=p(3);
y0=p(4);
sigx=p(5);
sigy=p(6);

z = b+ a*exp(-0.5 *( ((x(:,1)-x0).^2)/(sigx^2)  + ((x(:,2)-y0).^2)/(sigy^2)));

