%%% sample script to call cluster_tetrode_cmd
memory
tank = 'adrec107_102610';
blocks = {'adrec107-22' 'adrec107-23' 'adrec107-24' 'adrec107-25' 'adrec107-26' 'adrec107-27' 'adrec107-28'}
select_ica = 'n';
subtract_mean = 'y';
n_ica = 8;
max_snips = 50000;
basename = '1-7';
output_path = 'D:\data\TDT\adrec107-analysis\Blocks22-28';
disp_figs = 'n';
cluster_tetrode_cmd(tank,blocks,select_ica,subtract_mean,n_ica,max_snips, basename,output_path, disp_figs);
memory
exit
