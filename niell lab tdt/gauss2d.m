function z = gauss2d(p,x);
a= p(1);
b=p(2);
x0=p(3);
y0=p(4);
sigx=p(5);
sigy=p(6);

z = b+ a*exp(-0.5 *( ((x(:,1)-x0).^2)/(sigx^2)  + ((x(:,2)-y0).^2)/(sigy^2)));

