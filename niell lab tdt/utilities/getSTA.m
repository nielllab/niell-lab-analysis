function sta = getSTA(clustfile,afile, block,redo)
%%% read in single unit spike times for a given block
%%% this is mostly just a matter of sorting the spikes from one block
%%% but nice to do it just in one line!


load(clustfile,'Block_Name','Tank_Name','frameEpocs');
blocknum = find(strcmp(block,Block_Name));
Block_Name = Block_Name{blocknum}

clear sta
load(afile,'sta');

movieFile = 'C:\Users\Angie Michaiel\Desktop\movie files\cortex\wn_cortex_012alpha1_5hzLg30Hz.mat';
     load(movieFile); 

if ~exist('sta','var') | length(sta)<blocknum   | redo

    
    
    wn = noise_analysis_angie(clustfile,afile,movieFile,block);
    [params all_fit all_img] = fit2dgabor_angie(wn);
   
    for i = 1:length(params);
       sta(blocknum).nx(i) = params(i).nx;
       sta(blocknum).ny(i) = params(i).ny;
       sta(blocknum).sigx(i) = params(i).sigx;
       sta(blocknum).sigy(i) = params(i).sigy;
       sta(blocknum).exp_var(i) = params(i).exp_var;
       sta(blocknum).all_fit = all_fit;
       sta(blocknum).params = params;
       sta(blocknum).all_img = all_img;
  
    save(afile,'sta','-append')
    end
end

sta = sta(blocknum);
