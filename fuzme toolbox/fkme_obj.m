function fobj=fkme_obj(alfa,Uereq,nclass,data,phi,maxiter,distype,toldif)
% objective function of fkme, to find alfa that satisty mean of extragrade
% membership
if(alfa==0),
    fobj=-Uereq;
elseif(alfa==1),
	fobj=1-Uereq;
else
    ndata=size(data,1);
    Uinit= initmember(0.1,nclass,ndata);
    [U, Ue, centroid, dist, W, obj] = fkme(nclass,data,Uinit,phi,alfa,maxiter,distype,toldif);
    Ueav=mean(Ue);
    fobj=(Uereq-Ueav);
end;