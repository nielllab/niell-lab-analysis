%% resave the data, I want one session per directory to keep things organized
NRKO_resavedata

%% this is required by all functions
NRKO_global % always run at the start of analyses
NRKO_lev0_readdata_peristim

%% spike triggered spectrum 
NRKO_global                                 % create the info file
NRKO_lev1_sts_convol_peristim               % compute the spike triggered LFP spectra according to the convolution algorithm
NRKO_lev2_ppc_convol_peristim               % compute the PPC values
NRKO_lev3_ppc_convol_group_peristim         % pool the PPC values for group analysis
NRKO_lev4_ppc_convol_group_plot_peristim    % produce plots on the group PPC and get statistics in place

%% compute ERPs around stimulus onset
NRKO_global
NRKO_lev1_lfp_erp_peristim

%% compute TFRs of the LFP around stimulus onset
NRKO_global
NRKO_lev1_lfp_tfr_peristim
%%

