%%% script to run bars and drift again, to fix :
%%%  1) use spontaneous rate for baseline in bars
%%%  2) don't subtract minumum value of tuning curve in fit_tuningcurve for
%%%  in drifting gratings

clear all

barsweep16d_cluster_spont
%drift_orientfreq_new

[fname, pname] = uigetfile('*.mat','analysis');
oldpname = pname;
load(fullfile(pname,fname));


bars_amp;
theta_ind = 0:pi/8:15*pi/8;
 i=sqrt(-1);
 
 for cell_n=1:size(cells,1)
    mu = (sum(bars_amp(cell_n,:).*exp(i*theta_ind)))/sum(bars_amp(cell_n,:));
    if isnan(mu);
        mu=0;
    end
    bars_dsi(cell_n) = abs(mu);
end

clear alldata
alldata(:,1:2) = cells;
alldata(:,3) = bars_OSI;
alldata(:,4) = bars_dsi;
alldata(:,5) = driftorient_osi;
alldata(:,6) = driftorient_dsi;
xlswrite('circvar',alldata);
alldata

clear alldata
%%% bars 16d
%%% responsiveness
prefamp = (bars_B' + bars_A1');
%alldata(:,1) = prefamp./(prefamp + bars_spont);
%%% orientation selectivity
alldata(:,1) = (prefamp-bars_null')./(prefamp + bars_null');
%%% preferred orientation
alldata(:,2) = bars_theta'*180/pi;
%%% orientation widht
alldata(:,3) = abs(bars_w)*180/pi;
%%% direction selectivity
alldata(:,4) = (bars_A1'-bars_A2')./(bars_A1'+bars_A2' + 2*bars_B');
%%% peak response
alldata(:,5) = prefamp;
%%% spontaneous
alldata(:,6) = bars_spont;
alldata(:,7) = rf_width*30;

xlswrite('bars',alldata);
alldata

