function d = dtheta2(theta1,theta2);
d = theta1-theta2;

d = mod(d+pi,2*pi)-pi;

    