pathname = '\\niell-v2-w7\Angie_analysis\DOI_experiments';
n=0;

n=n+1;
files(n).expt = '120215';  %%% can be any useful name
files(n).dir = '12_02_15';  %%% directory within main path that has data for this expt
files(n).tank = 'tankname';  %%% enter name of the tank (as registered in TDT)
files(n).clusterfile = '120215cluster';
files(n).analysisfile = '120215analysis';
files(n).darkness = 'dark';   %%% name of block for different stim types
files(n).gratings = 'drift'
files(n).wn = 'noise';
files(n).bars = 'bars1';
files(n).movespots = '';  %%% leave empty if there is no data
files(n).treatment = 'DOI';
files(n).layers = [5 6]; 
files(n).notes = 'good data';
files(n).misc = 'nice units! but some movement artifact';

n=n+1;
files(n).expt = '091815';
files(n).dir = '09_18_15';
files(n).tank = 'tankname';  %%% enter name of the tank (
files(n).clusterfile = '';
files(n).analysisfile = '';
files(n).darkness = 'dark';   %%% name of block for different stim types
files(n).gratings = 'drift'
files(n).wn = 'noise';
files(n).bars = 'bars1';
files(n).movespots = '';  %%% leave empty if there is no data
files(n).treatment = 'DOI';
files(n).layers = [5 6]; 
files(n).notes = 'good data';
files(n).misc = 'nice units! but some movement artifact';