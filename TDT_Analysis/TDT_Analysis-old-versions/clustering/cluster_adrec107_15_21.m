%%% sample script to call cluster_tetrode_cmd
memory
tank = 'adrec107_102610';
blocks = {'adrec107-15' 'adrec107-16' 'adrec107-17' 'adrec107-18' 'adrec107-19' 'adrec107-20' 'adrec107-21'}
select_ica = 'n';
subtract_mean = 'y';
n_ica = 8;
max_snips = 50000;
basename = '1-7';
output_path = 'D:\data\TDT\adrec107-analysis\Blocks15-21';
disp_figs = 'n';
cluster_tetrode_cmd(tank,blocks,select_ica,subtract_mean,n_ica,max_snips, basename,output_path, disp_figs);
memory
exit
