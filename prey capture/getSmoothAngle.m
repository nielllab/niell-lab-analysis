function theta = getSmoothAngle(vec)
headtheta = atan2(vec(:,1),vec(:,2));
vtheta = diff(headtheta)
vtheta(vtheta>pi) = vtheta(vtheta>pi)-2*pi;
vtheta(vtheta<-pi) = vtheta(vtheta<-pi)+2*pi;
theta = [0 ; cumsum(vtheta)] + headtheta(1);