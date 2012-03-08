%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sample script to call cluster_tetrode_cmd

clear all
pack

tank = '022309_awake_tet';
blocks = {'bars16d1' 'drift1'}
select_ica = 'n';
subtract_mean = 'y';
n_ica = 8;
max_snips = 50000;
basename = 'bardrift1';
output_path = 'C:\Documents and Settings\idl\My Documents\MATLAB\0223_testnofigs';
disp_figs = 'n';
cluster_tetrode_cmd(tank,blocks,select_ica,subtract_mean,n_ica,max_snips, basename,output_path, disp_figs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
