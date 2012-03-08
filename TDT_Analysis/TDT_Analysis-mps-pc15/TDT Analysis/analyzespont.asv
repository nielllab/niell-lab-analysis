function analyzespont
% Matlab codes for reading EEG data
% For each condition, goes by epoc and plots power spectrum
% cmn 06-06, based on code by Jianhua Cang 06-27-03



%%% read in cluster data, (only needed to get block name and start/stop time)
%%% then connect to the tank and read the block

%%% read in cluster data, then connect to the tank and read the block
clear all
cells =1;
[fname, pname] = uigetfile('*.mat','cluster data');
load(fullfile(pname,fname));
for i =1:length(Block_Name);
    sprintf('%d : %s ',i,Block_Name{i})
end
block = input('which block to analyze ? ');
Block_Name = Block_Name{block}

spontnum = input('which spont to store as ? ')


TTX = openTTX(Tank_Name,Block_Name); % to initialize
time1= 0;
time2=0;
max_events = 50000;

[tsamp vsmooth] = getBlockVelocity(Tank_Name,Block_Name);

thresh_velocity=2;

[afname, apname] = uigetfile('*.mat','analysis data');
if afname~=0
    afile_spont = fullfile(apname,afname);
    save(afile_spont, 'afile_spont','apname','-append'); %%% to make sure it's not overwritten
    load(afile_spont);
    use_afile=1;
else
    use_afile=0;
    %%% select which units to analyze (channel, cluster number)

end


[tsamp vsmooth groomstate] =getBlockVelocity_groom(Tank_Name,Block_Name);

dt = mean(diff(tsamp))
 
    for c = 1:size(cells,1)


        ch = cells(c,1);
        clust = cells(c,2);

block_times = event_times_all-(block-1)*10^5;
        channel_times =squeeze(block_times(ch,:));
        used = find(idx_all(ch,:)==clust & channel_times>0 & channel_times<10^5);

        spiketimes = channel_times(used);
        spiketimes = spiketimes(spiketimes>min(tsamp) & spiketimes<max(tsamp));

        
        v_interp = interp1(tsamp,vsmooth,spiketimes);
        g_interp = interp1(tsamp,groomstate,spiketimes);
        r(c,1) = sum(v_interp<1)/(dt*sum(vsmooth<1));
        r(c,2) = sum(v_interp>2)/(dt*sum(vsmooth>2));
    
        spontrate{spontnum}=r;
        spontrate

end %% tet
figure
bar(r)
save(afile_spont, 'spontrate' ,'-append');
invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');