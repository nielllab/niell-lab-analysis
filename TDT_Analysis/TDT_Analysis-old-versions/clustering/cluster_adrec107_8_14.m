%%% sample script to call cluster_tetrode_cmd
memory
tank = 'adrec107_102610';
blocks = {'adrec107-8' 'adrec107-9' 'adrec107-10' 'adrec107-11' 'adrec107-12' 'adrec107-13' 'adrec107-14'}
select_ica = 'n';
subtract_mean = 'y';
n_ica = 8;
max_snips = 50000;
basename = '1-7';
output_path = 'D:\data\TDT\adrec107-analysis\Blocks8-14';
disp_figs = 'n';
cluster_tetrode_cmd(tank,blocks,select_ica,subtract_mean,n_ica,max_snips, basename,output_path, disp_figs);
memory
exit
