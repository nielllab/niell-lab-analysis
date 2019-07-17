function newAngles = wrapAng(oldAng)

dangle = diff(oldAng);
dangle(dangle>180) = dangle(dangle>180)-360;
dangle(dangle<-180) = dangle(dangle<-180)+360;
newAngles = oldAng(1) + cumsum(dangle,'omitnan');
newAngles = [oldAng(1); newAngles];