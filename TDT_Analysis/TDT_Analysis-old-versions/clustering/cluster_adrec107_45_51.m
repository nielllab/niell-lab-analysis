%%% sample script to call cluster_tetrode_cmd
memory
tank = 'adrec107_102610';
blocks = {'adrec107-45' 'adrec107-46' 'adrec107-47' 'adrec107-48' 'adrec107-49' 'adrec107-50' 'adrec107-51' }
select_ica = 'n';
subtract_mean = 'y';
n_ica = 8;
max_snips = 50000;
basename = '1-7';
output_path = 'D:\data\TDT\adrec107-analysis\Blocks45-51';
disp_figs = 'n';
cluster_tetrode_cmd(tank,blocks,select_ica,subtract_mean,n_ica,max_snips, basename,output_path, disp_figs);
memory
exit
