%%appending Indices or variables that you want to analyze at the end i.e inh waveforms versus exc, there may be a cleaner
%%way to do this
append_analysis_file

%% extracting the data from Niell lab analysis and cluster files
getBlockData_specific_block_sessions

%% resave the data, I want one session per directory to keep things organized
NRKO_resavedata

%% this is required by all functions
NRKO_global % always run at the start of analyses will be specific to project
NRKO_lev0_readdata_peristim %converts to field trip format and filters out line noise from LFP spec

%% produce peristimulus spike times
NRKO_global
NRKO_lev1_sdf_peristim%last run 3-19-16 on 1-26-16 data xtraction 
% added in tuning as calculated by Cris's programs instead, mesh with
% martin's analysis below?
NRKO_lev2_tuningstats_peristim
NRKO_lev3_tuningstats_group_peristim
NRKO_lev4_tuning_group_plot_peristim

%% spike triggered spectrum 
NRKO_global                                 % create the info file
NRKO_lev1_sts_convol_peristim               % compute the spike triggered LFP spectra according to the convolution algorithm
NRKO_lev2_ppc_convol_peristim               % compute the PPC values
%NRKO_lev2_ppc_convol_peristim_CSD           % compute the PPC values relative to the CSD, this didn't change much in the locking measures in the NRKO mice 
%NRKO_lev3_ppc_convol_group_peristim_CSD
NRKO_lev3_ppc_convol_group_peristim         % pool the PPC values for group analysis
NRKO_lev4_ppc_convol_group_plot_peristim    % produce plots on the group PPC and get statistics in place
%NRKO_lev4_ppc_convol_group_plot_peristim_CSD
NRKO_lev4_spec_group_plot_peristim

%% 
NRKO_global
NRKO_lev1_lfp_spectrum_peristim %last run 1-26-16
NRKO_lev2_lfp_spectrum_granger_peristim_group
NRKO_lev2_lfp_spectrum_peristim_group
% NRKO_lev2_spec_peristim
%% compute ERPs around stimulus onset
% NRKO_global
% NRKO_lev1_lfp_erp_peristim
% 
% %% compute TFRs of the LFP around stimulus onset
% NRKO_global
% NRKO_lev1_lfp_tfr_peristim
% %%

