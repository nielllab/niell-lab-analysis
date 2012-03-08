clear all
cells=0;
afiles = {'C:\data\matlab data\data 2009\021609_awake\allstim1_new\analysis.mat' ...
    'C:\data\matlab data\data 2009\021609_awake\allstim2_new\analysis2.mat' ...
    'C:\data\matlab data\data 2009\030609_awake_tet_analysis\allstim1\analysis.mat', ...
    'C:\data\matlab data\data 2009\022209_awake_linear_analysis\allstim3\analysis.mat', ....
    'C:\data\matlab data\data 2009\022309_awake_tet_analysis\allstim1\analysis.mat' ...
    'C:\data\matlab data\data 2009\022309_awake_tet_analysis\allstim2\analysis.mat'};
N =0;
for i = 1:length(afiles)
 
    load(afiles{i});
    n_units = length(L_ratio);
    cellrange = N+1:N+n_units;
    N=N+n_units;
    
    driftA1(cellrange,:) =driftorient_A1;
    driftA2(cellrange,:) = driftorient_A2;
    driftB(cellrange,:) = driftorient_B;
    drift_theta_w(cellrange,:)=driftorient_thetawidth;
    driftspont(cellrange,:) = drift_spont;
    
end
