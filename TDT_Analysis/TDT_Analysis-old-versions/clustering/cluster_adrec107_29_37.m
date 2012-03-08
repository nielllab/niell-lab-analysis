%%% sample script to call cluster_tetrode_cmd
memory
tank = 'adrec107_102610';
blocks = {'adrec107-29' 'adrec107-30' 'adrec107-31' 'adrec107-32' 'adrec107-33' 'adrec107-34' 'adrec107-35'  'adrec107-36'  'adrec107-37'}
select_ica = 'n';
subtract_mean = 'y';
n_ica = 8;
max_snips = 50000;
basename = '1-7';
output_path = 'D:\data\TDT\adrec107-analysis\Blocks29-37';
disp_figs = 'n';
cluster_tetrode_cmd(tank,blocks,select_ica,subtract_mean,n_ica,max_snips, basename,output_path, disp_figs);
memory
exit
