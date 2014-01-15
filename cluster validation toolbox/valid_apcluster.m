% Affinity Propagation clustering
nrun = 2000;     % max iteration times, default 2000
nconv = 50;      % convergence condition, default 50
lam = 0.90;      % damping factor, default 0.9
splot = 'noplot';
%splot = 'plot'; % observing a clustering process when it is on

M = simatrix_make(data,dtype,nrow);

k1 = round(2*sqrt(nrow))-1;
[labelid,netsim,dpsim,expref,pref,lowp,highp]=apclusterK(...
    M,1,nrun,nconv,lam,splot,0,0,0,1,k1);

% AP gives k clusters exactly under given k
[labelid,netsim,dpsim,expref,pref]=apclusterK(...
    M,k,0,nrun,nconv,lam,splot,lowp,highp,1,k1);

% NC = length(unique(labelid));
[C, label] = ind2cluster(labelid);
