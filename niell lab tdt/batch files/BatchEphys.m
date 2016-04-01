% file entry info for batch processing of single unit analysis E-phys
% data with or without optogenetic tagging

%pathname = 'D:\Jen_analysis\NR5A_Pinping\';

n=0;

%%%% 
% n=n+1;
% 
% % subject parameters
% 
% files(n).condition=1;%control =1;mutant=2; 
% files(n).sex=1; %1=female and 2=male
% files(n).age=1;%1=adult, 2=p21, 3= EO3, 4=EO1 etc.
% 
% %file parameters and paths
% files(n).expt = '';
% files(n).path = 'D:\Jen_analysis\NR5A_Pinping\6_17_15\';
% files(n).clusterfile = 'cluster_data_6_17_15_6_17_15.mat';
% files(n).analysisfile = 'analysis_2.mat';
% files(n).pathtank='D:\Jen_tanks\6_17_15\'
% files(n).saveFile='drift_batch';
% files(n).notes = '';
% 
% %block info
% files(n).blockDrift = {'drift1_1secISI'};
% files(n).blockWn = {'wn1'};
% files(n).blockPinp = {'pimp1'};
% files(n).blockBar = {'bars1' };
% files(n).postBlocks={''};
% files(n).nchan=64; %number of sites
% 
% %parameters for layer analysis
% files(n).tip_loc_1=550;
% files(n).tip_loc_2=550;
% files(n).angle=45; 

%%
n=n+1;
files(n).condition=1;%control =1;mutant/drug/etc.=2; 
files(n).sex=1; %1=female and 2=male
files(n).age=1;%1=adult, 2=p21, 3= EO3, 4=EO1 etc.
files(n).expt = '';%could enter date in order to organize by date
files(n).notes = ''; %could comment whether experiment was good bad or other flags

%file parameters and paths

files(n).path = 'D:\Angie_analysis\DOI_experiments\01_08_16\';
files(n).clusterfile = 'cluster_data_01_08_16htr2a_doi_cluster_010816.mat';
files(n).analysisfile = 'analysis_010816.mat';
files(n).pathtank='D:\Angie_tanks\01_08_16htr2a_doi\';

%block info
files(n).blockDrift = {'drift_predoi1' 'drift_postdoi1'};
files(n).blockWn = {'wn_predoi1' 'wn_postdoi1'};
files(n).blockPinp = {'laser_predoi_int05' 'laser_pre_int_075' 'laser_postdoi_int075' };
files(n).blockBar = {''};
files(n).blockNoStim={''};
files(n).postBlocks={'dark_postdoi1' 'drift_postdoi1' 'laser_postdoi_int075' 'wn_postdoi1'};
files(n).nchan=64; %number of sites

%parameters for layer analysis
files(n).tip_loc_1=550;
files(n).tip_loc_2=550;
files(n).angle=45; 







