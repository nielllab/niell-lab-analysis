

clear all

dbstop if error
n=0;

%%% 'analysis_11302012_rec1c.mat'...
%%% 'analysis_03122013_rec1.mat'...

   
analysisPath = 'C:\data\matlab data\Wayne Matlab data\analysis files\awake lgn\'
afile = {     'analysis_11282012_rec1c.mat'...
    'analysis_11282012_rec2c.mat'...
    'analysis_11282012_rec3c.mat'...
    'analysis_11292012_rec1c.mat'...
    'analysis_11292012_rec3c.mat'...
   'analysis_11302012_rec2c.mat'...
    'analysis_12062012_rec1c.mat'...
    'analysis_12072012_rec1b.mat'...
    'analysis_12072012_rec2c.mat'...
    'analysis_12072012_rec3c.mat'...
    'analysis_01092013_rec3.mat'...
    'analysis_01162013_rec3.mat'...
    'analysis_01232013_rec3.mat'...
    'analysis_01242013_rec2.mat'...
    'analysis_02182013_rec1.mat'...
    'analysis_02182013_rec2.mat'...
    'analysis_02192013_rec1.mat'...
    'analysis_02192013_rec2.mat'...
    'analysis_02202013_rec1.mat'...
    'analysis_02212013_rec1.mat'...
    'analysis_02212013_rec2.mat'...
    'analysis_0307a2013_rec1.mat'...
    'analysis_0307a2013_rec2.mat'...
    'analysis_03082013_rec1.mat'...
    'analysis_03112013_rec1.mat'...
    'analysis_03122013_rec2.mat'...
    };


extraPath = 'C:\data\matlab data\Wayne Matlab data\analysis files\awake lgn new\'
for fnum = 1:length(afile)
    %for fnum=26:26
    close all
    fnum
    fullfile(extraPath,[afile{fnum}(1:end-4) '_new.mat'])
    fullfile(extraPath,[afile{fnum}(1:end-5) '_new.mat'])
    
   if exist(fullfile(extraPath,[afile{fnum}(1:end-4)  '_new.mat']),'file') | exist(fullfile(extraPath,[afile{fnum}(1:end-5) '_new.mat']),'file')
     if exist(fullfile(extraPath,[afile{fnum}(1:end-4) '_new.mat']),'file')   
       load(fullfile(extraPath,[afile{fnum}(1:end-4) '_new.mat']))
   elseif exist(fullfile(extraPath,[afile{fnum}(1:end-5) '_new.mat']),'file')
     load(fullfile(extraPath,[afile{fnum}(1:end-5) '_new.mat']))
     end
     
     useCells=cells

        
        analysisFile = fullfile(analysisPath,afile{fnum})
        load(analysisFile);
         cells
        clusterfilename
        pname
        
        if exist(clusterfilename,'file')
            clusterFile = clusterfilename;
        elseif exist(clusterfilename((length(pname)+1):end),'file')
            clusterFile = clusterfilename((length(pname)+1):end);
        elseif exist([clusterfilename((length(pname)+1):end) '.mat'],'file')
            clusterFile = [clusterfilename((length(pname)+1):end) '.mat'];
        elseif exist([clusterfilename '.mat'],'file')
            clusterfile = [clusterfilename '.mat'];
            
        else
           keyboard
           [fname pname] = uigetfile('*.mat','cluster file');
            clusterFile = fullfile(pname,fname);
            clusterfilename = clusterFile;
            if fname~=0
                save(analysisFile,'clusterfilename','-append');
            end
        end
        clusterFile
        load(clusterFile);
        blocknum=0;
       
        
        for blocknum = 1:length(Block_Name)
            
            
            blocknum
            Block_Name{blocknum};
            
            flags = struct('visStim',1,'lfpTseries',1,'lfpSpectra',1,'mouseOn',1);
            % flags = struct('visStim',1);
            
            ch = unique(useCells(:,1))';
            
            data = getTDTdata(Tank_Name,Block_Name{blocknum},ch,flags);
            if isnan(data.frameEpocs)
                sprintf('check if tank %s is registered',Tank_Name)
                break
            end
            exptdata{fnum}.mouseT{blocknum} = data.mouseT;
            exptdata{fnum}.mouseV{blocknum} = data.mouseV;
            exptdata{fnum}.block{blocknum} = Block_Name{blocknum};
            exptdata{fnum}.tank = Tank_Name;
            exptdata{fnum}.analysis_file = afile{fnum};
            exptdata{fnum}.cluster_file = clusterFile;
            exptdata{fnum}.stimEpocs{blocknum}=data.stimEpocs;
            exptdata{fnum}.frameEpocs{blocknum} = data.frameEpocs;
            for i = 1:length(ch)
                exptdata{fnum}.lfpT{blocknum}{i} = data.lfpT{ch(i)};
            exptdata{fnum}.lfpV{blocknum}{i} = data.lfpData{ch(i)};
            end
            exptdata{fnum}.lfpChans = ch;
            %         exptdata{fnum}.displayHeight = displayHeight;
            %         exptdata{fnum}.displayOffset = displayOffset;
            for cell_n = 1:size(useCells,1);
                matchCell = find(cells(:,1) == useCells(cell_n,1) & cells(:,2)==useCells(cell_n,2));
                if isempty(matchCell)
                    sprintf('couldnt find match ch=%d cl=%d',useCells(:,1),useCells(:,2))
                else
                    unitdata{cell_n+n}.spikes{blocknum} =spikeT{matchCell}(spikeT{matchCell}>(blocknum-1)*10^5 & spikeT{matchCell}<(blocknum-1 + 0.5)*10^5) - (blocknum-1)*10^5;
                    unitdata{cell_n+n}.expnum = fnum;
                    s= unitdata{cell_n+n}.spikes{blocknum};
                    tic
                  %  [unitdata{cell_n}.coherence{blocknum} t] = spikeLFPcoh(s,data.lfpT{useCells(cell_n,1)},data.lfpData{useCells(cell_n,1)},0.004,0.5);
                   % [unitdata{cell_n}.coh_shuffle{blocknum} t] = spikeLFPcoh(rand(length(s),1)*max(s),data.lfpT{useCells(cell_n,1)},data.lfpData{useCells(cell_n,1)},0.004,0.5);
                    unitdata{cell_n}.ch = useCells(cell_n,1);
                    unitdata{cell_n}.clust = useCells(cell_n,2);
                    
%                     d = diff(s);
%                     figure
%                     hist(d(d<0.05),0:0.001:0.05);
%                     toc
%                     figure
%                     plot(t,unitdata{cell_n}.coherence{blocknum})
%                     hold on
%                     plot(t,unitdata{cell_n}.coh_shuffle{blocknum},'g')
                end
                
            end
        end
        n=n+size(useCells,1);
   else
        sprintf('couldnt find file %d %s',fnum,afile{fnum})
    end
end

% manual_type = xlsread('C:\data\lgn rf project\lgn_analysis\lgn types 091412.xlsx','A1:A294');
% manual_outside = xlsread('C:\data\lgn rf project\lgn_analysis\lgn types 091412.xlsx','B1:B294');
% manual_outside(isnan(manual_outside))=0;
%
% used =(manual_type~=0);

% load 'C:\data\lgn rf project\lgn_analysis\LGNunits'
% use_list=find(used);
%  for i = 1:length(use_list);
%      unitdata_select{i} = unitdata{use_list(i)};
%  end
%  unitdata=unitdata_select;

%
save('Tschetter_adult_awake032614','exptdata','unitdata','-v7.3');