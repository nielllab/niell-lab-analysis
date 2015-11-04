

clear all

dbstop if error
n=0;

   
analysisPath = 'D:\Jen_analysis\NR5A_Pinping\Jen_NR5A_analysis_files\analysis_files\'

for dataset =1:3
   
    if dataset==1
    afile = {'5_4_15_analysis_2.mat',...
    '5_11_15_analysis_2.mat',...
   '5_14_15_analysis_2.mat',...
   '6_25_15_analysis_2.mat',...
   '6_29_15_analysis_2A.mat',...
   '8_10_15_analysis_2.mat',...
   '8_17_15_analysis_2.mat'}
%     'analysis_12_22_14_1st_clust.mat'}; no drift session
%'2_25_15_analysis_2.mat',... shorter ISI during drift session
  %   '3_11_15_analysis_2.mat',...
%   '4_9_15_analysis_2.mat',...
%     '4_10_15_analysis_2.mat',...
%'1_13_15_analysis_pos_ctl.mat'
%'4_13_15_analysis_2.mat',...
%'4_30_15_analysis_2.mat',...
   
    elseif dataset==2
        
    analysisPath = 'D:\Jen_analysis\NR5A_Pinping\Jen_NR5A_analysis_files\analysis_files\N2B\'

    afile = { '4_23_15_analysis_2.mat',...
    '6_18_15_analysis_2A.mat',...
    '6_22_15_analysis_2A.mat',...
    '6_23_15_analysis_2.mat',...
    '6_26_15_analysis_2A.mat',...
    '7_1_15_analysis_2.mat',...
    '8_11_15_analysis_2A.mat',...
    '8_11_15_rec2_analysis_2.mat',...
    '8_13_15_analysis_2A.mat',...
    '8_13_15_rec2_analysis_2A.mat',...
    '8_14_15_analysis.mat',...
    '8_19_15_analysis_2.mat',...
    '8_19_15_rec2_analysis_2.mat'};
%'2_24_15_analysis_2.mat',...bars corrupt?
%'3_25_15_analysis_3_25_15.mat',...
%  '3_26_15_analysis_2.mat',...
 %   '2_2_15_analysis_2.mat',...
    else dataset==3
        analysisPath = 'D:\Jen_analysis\NR5A_Pinping\Jen_NR5A_analysis_files\analysis_files\N2A\'

    afile = {'5_12_15_analysis_2.mat',...
     '6_28_15_analysis_2A.mat',...
   '6_30_15_analysis_2B.mat',...
    '7_29_15_analysis_2.mat',...
   '8_7_15_ analysis_2A.mat',...
   '8_18_15_analysis_2.mat',...
    '8_18_15_rec2_analysis_2.mat',...
    '9_12_15_analysis.mat',...
    'analysis_9_14_15.mat',...
    '9_16_15_analysis_2.mat'};
%'8_6_15_analysis_2.mat',...bar session corrupted
%'3_3_15_analysis_2.mat'...
      %  '3_4_15_analysis_2.mat',...
     %   '4_24_15_analysis_2A.mat',...
    end



for fnum = 1:length(afile)
    %for fnum=26:26
    
    
    close all
    clear B
    
        fnum        
        analysisFile = fullfile(analysisPath,afile{fnum})
        load(analysisFile);
         
        cells
        
%         clear  clusterfilename %replace clusterfile names generated on   other computers
        pname
        useCells=cells
        afile(fnum)

        if exist(clusterfilename,'file') 
             clusterFile = clusterfilename;    
        else
           %keyboard
           [fname pname] = uigetfile('*.mat','cluster file');
            clusterFile = fullfile(pname,fname);
            clusterfilename = clusterFile;
            if fname~=0
                save(analysisFile,'clusterfilename','-append');
            end
        end
            
        afile(fnum)
        clusterFile
        load(clusterFile);
         
%         for blocknum = 3:4
            Block_Name
            B = strncmpi(Block_Name,'d',1) 
            
                 if sum(B)>1 % && B(1,2)~=1   %or not equal to the whatever the most commone stim order is
            blocknum = listdlg('ListString',Block_Name,'SelectionMode','single');
                 else
            blocknum=find(B==1);
                 end
                  
            if sum(B)>=1
                      
            Block_Name{blocknum};
            
            flags = struct('visStim',1,'lfpTseries',1,'lfpSpectra',1,'mouseOn',1);
            % flags = struct('visStim',1);
            
            ch = unique(useCells(:,1))';
          
            
            data = getTDTdata(Tank_Name,Block_Name{blocknum},ch,flags);
            if isnan(data.frameEpocs)
                sprintf('check if tank %s is registered',Tank_Name)
%                 return
            end
  
            exptdata{fnum,dataset}.mouseT{1} = data.mouseT;
            exptdata{fnum,dataset}.mouseV{1} = data.mouseV;
            exptdata{fnum,dataset}.block{1} = Block_Name{blocknum};
            exptdata{fnum,dataset}.tank = Tank_Name;
            exptdata{fnum,dataset}.analysis_file = afile{fnum};
            exptdata{fnum,dataset}.cluster_file = clusterFile;
            exptdata{fnum,dataset}.stimEpocs{1}=data.stimEpocs;
            exptdata{fnum,dataset}.frameEpocs{1} = data.frameEpocs;
            
           % keyboard
            
            for i = 1:length(ch)
                exptdata{fnum,dataset}.lfpT{1}{i} = data.lfpT{ch(i)};
            exptdata{fnum,dataset}.lfpV{1}{i} = data.lfpData{ch(i)};
            end
            
            exptdata{fnum,dataset}.lfpChans = ch;
            
            
                for cell_n = 1:size(useCells,1);
               
                    matchCell = find(cells(:,1) == useCells(cell_n,1) & cells(:,2)==useCells(cell_n,2));
                if isempty(matchCell)
                    sprintf('couldnt find match ch=%d cl=%d',useCells(:,1),useCells(:,2))
                else
                    
                    unitdata{cell_n+n}.spikes{1} =spikeT{matchCell}(spikeT{matchCell}>(blocknum-1)*10^5 & spikeT{matchCell}<(blocknum-1 + 0.5)*10^5) - (blocknum-1)*10^5;
                    unitdata{cell_n+n}.expnum = fnum;
                    unitdata{cell_n+n}.GT = dataset;
                    unitdata{cell_n+n}.block = Block_Name{blocknum};
                    unitdata{cell_n+n}.pinp = pinped(cell_n);
                    unitdata{cell_n+n}.responsive = resp(cell_n);
                    unitdata{cell_n+n}.layer=layer(cell_n);
                    unitdata{cell_n+n}.inhibitory = inh(cell_n);
                    s= unitdata{cell_n+n}.spikes{1};
                    tic
                  
                    unitdata{cell_n +n}.ch = useCells(cell_n,1);
                    unitdata{cell_n +n}.clust = useCells(cell_n,2);
        
                end
                
                end
%         end
                    n=n+size(useCells,1);
                  end
                  
end
  % else
      %  sprintf('couldnt find file %d %s',fnum,afile{fnum})
   % end
end


save('JLH_NMDA_KO_drift_test1','exptdata','unitdata','-v7.3');