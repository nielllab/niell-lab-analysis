function d = dtheta_old(theta1,theta2);
d = theta1-theta2;

while min(d)<-pi
    d(d<-pi)=d(d<-pi)+2*pi;
end
while max(d)>pi
    d(d>pi) = d(d>pi)-2*pi;
end

    