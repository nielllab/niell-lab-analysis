%F:\Jennifer_Development\Initial_development_tanks\Jen_developmental_Tanks

dbstop if error
n=0;

   
analysisPath = 'F:\Jennifer_Development\';

for dataset =1:2
   
    if dataset==1
    afile = {'Good recordings\Adults\5_25_13_Adultmale3mo\Adult B\analysis_5_25_13_B_rec1_adult.mat',...
            'Good recordings\Adults\7_25_13_deep\MouseB\Rec1\analysis_7_25_13_mouseB_rec1.mat',...
            'Good recordings\Adults\7_25_13_deep\MouseB\Rec2\analysis_cluster_data_07_25_13_mouseB_adult_rec2.mat',...
            'Good recordings\Adults\9_12_13\rec1\analysis_9_12_13_adult_rec1.mat',...
            'Good recordings\Adults\11_11_13\rec1\analysis_11_11_13_adult_rec1.mat',... %%% no wn
            'Good recordings\Adults\11_11_13\rec2\analysis_11_11_13_adult_rec2.mat',...
            'Good recordings\Adults\11_13_13\rec1\analysis_adult_11_13_13_rec1.mat',...
            'Good recordings\Adults\11_13_13\rec2\analysis_11_13_13_rec2.mat',...
            'Good recordings\Adults\11_14_13\rec1\analysis_11_14_13_adult_rec1.mat',...
            'Good recordings\Adults\11_14_13\rec2\analysis_11_14_13_adult_rec2.mat',...
             'Good recordings\Adults\11_15_13\rec1\analysis_11_15_13_adult_rec1.mat',...
             'Good recordings\Adults\11_15_13\rec2\analysis_11_15_13_adult_rec2.mat'}
 %'Good recordings\Adults\7_19_13_adult\analysis_cluster_adult_rec2_7_19_13.mat',...
   %'Good recordings\Adults\1month_2month old\4_22_13\analysis_adult.mat',...
   %'Good recordings\Adults\1month_2month old\9_25_13\rec1\analysis_9_25_13_P31_rec1.mat',...
    elseif dataset==2
        
    afile = {'Good recordings\EO1_EO2\8_7_13_EO1\rec1_full_clustering\analysis_8_7_13_EO1_rec1_more_strigent.mat',...
            'Good recordings\EO1_EO2\8_7_13_EO1\rec2_full_clustering\analysis_8_7_13_rec2.mat',... 
            'Good recordings\EO1_EO2\8_8_13_EO2\rec1_full_clustering\analysis.mat',...
            'Good recordings\EO1_EO2\5_22_13_EO1\analysis_rec1_A_5_22_13_strict_selection.mat',...         
            'Good recordings\EO1_EO2\9_9_13_EO1\rec1\analysis_9_9_13_EO1_rec1.mat',...
            'Good recordings\EO1_EO2\9_9_13_EO1\rec2\analysis_9_9_13_EO1_rec2.mat',...
            'Good recordings\EO1_EO2\9_30_13_EO1\rec1\analysis_9_30_13_EO1_rec1.mat',...
            'Good recordings\EO1_EO2\9_30_13_EO1\rec2\analysis_9_30_13_EO1_rec2.mat',...
            'Good recordings\EO1_EO2\11_20_13_EO0\rec1\analysis_11_20_13_rec1.mat',...
            'Good recordings\EO1_EO2\11_20_13_EO0\rec2\analysis_11_20_13_EO0_rec2.mat',...
            'Good recordings\EO1_EO2\11_21_13_EO1_mouseA\analysis_11_21_13_mouseA_EO1.mat',...
            'Good recordings\EO1_EO2\11_23_13_EO2\analysis_11_23_13_EO2.mat',...
            'Good recordings\EO1_EO2\11_26_13_mouseA_EO1\rec1\analysis_11_26_13_EO1_mouseA_rec1.mat',...
            'Good recordings\EO1_EO2\11_26_13_mouseA_EO1\rec2\analysis_11_26_13_EO1_rec2.mat',...
            'Good recordings\EO1_EO2\11_26_13_mouseB_EO1\analysis_11_26_13_mouseB_EO1.mat'};
%'Good recordings\EO1_EO2\10_1_13_EO2\rec1\analysis_10_1_13_EO2_rec1.mat',...
%            'Good recordings\EO1_EO2\4_29_13_EO1\mouseC\analysis_4_29_13_C.mat',... 
%'Good
%recordings\EO1_EO2\7_17_13_EO1\Analysis_7_17_13_cluster_7_17_13_EO1.mat',...%bars
%not clustered

%     else dataset==3
%         analysisPath = 'D:\Jen_analysis\NR5A_Pinping\Jen_NR5A_analysis_files\analysis_files\N2A\'
% 
%     afile = {}   
    end



for fnum = 1:length(afile)
    %for fnum=26:26
    
    
    close all
    clear B
    
        fnum        
        analysisFile = fullfile(analysisPath,afile{fnum});
        load(analysisFile);
         
        cells
        
%         clear  clusterfilename %replace clusterfile names generated on   other computers
        pname
        useCells=cells;
        afile(fnum);

        if exist(clusterfilename,'file') ;
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
        clusterFile;
        load(clusterFile);
         
        clusterFile
%         for blocknum = 3:4
            Block_Name
            %B = strncmpi(Block_Name,'Bar',3) ;
            
%                  if sum(B)>1 ;
%             blocknum = listdlg('ListString',Block_Name,'SelectionMode','single');
%                  else
%             blocknum=find(B==1);
%                  end
%                   
 blocknum = listdlg('ListString',Block_Name,'SelectionMode','single');
 
           % if sum(B)~=0
                      
            Block_Name{blocknum};
            
            flags = struct('visStim',1,'lfpTseries',1,'lfpSpectra',1,'mouseOn',1);
            % flags = struct('visStim',1);
            
            ch = unique(useCells(:,1))';
          
            ch_lfp=1:64;
            tank_path='F:\Jennifer_Development\Initial_development_tanks\Jen_developmental_Tanks\'
            Tank_Name=fullfile(tank_path,Tank_Name) 
            
            data = getTDTdata(Tank_Name,Block_Name{blocknum},ch_lfp,flags);
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
            g=(useCells(:,1));
            depth=[g layer];
            exptdata{fnum,dataset}.layer{1}=depth;
           % keyboard
            
         
           
            for i = 1:64    
            exptdata{fnum,dataset}.lfpT{i} = data.lfpT{i};
            exptdata{fnum,dataset}.lfpV{i} = data.lfpData{i};    
            end
            
            exptdata{fnum,dataset}.lfpChans = i;
            
            
                for cell_n = 1:size(useCells,1);
               
                    matchCell = find(cells(:,1) == useCells(cell_n,1) & cells(:,2)==useCells(cell_n,2));
                if isempty(matchCell)
                    sprintf('couldnt find match ch=%d cl=%d',useCells(:,1),useCells(:,2))
                else
                    
                    unitdata{cell_n+n}.spikes{1} =spikeT{matchCell}(spikeT{matchCell}>(blocknum-1)*10^5 & spikeT{matchCell}<(blocknum-1 + 0.5)*10^5) - (blocknum-1)*10^5;
                    unitdata{cell_n+n}.expnum = fnum;
                    unitdata{cell_n+n}.GT = dataset;
                    unitdata{cell_n+n}.block = Block_Name{blocknum};
                    %unitdata{cell_n+n}.pinp = pinped(cell_n);
                    unitdata{cell_n+n}.responsive = resp(cell_n);
                    unitdata{cell_n+n}.layer=layer(cell_n);
                    unitdata{cell_n+n}.inhibitory = inh(cell_n);
                    unitdata{cell_n+n}.OSI = OS(cell_n);
                    %unitdata{cell_n+n}.cvOSI = cvOS(cell_n);
                    unitdata{cell_n+n}.DSI = DS(cell_n);
                    %unitdata{cell_n+n}.cvDSI = cvDS(cell_n);
                    unitdata{cell_n+n}.SF = SF(cell_n);
                    unitdata{cell_n+n}.Opref = prefOrient(cell_n);
                    unitdata{cell_n+n}.width = TW(cell_n);
                    s= unitdata{cell_n+n}.spikes{1};
                    tic
                  
                    unitdata{cell_n +n}.ch = useCells(cell_n,1);
                    unitdata{cell_n +n}.clust = useCells(cell_n,2);
        
                end
                
                end
%         end
                    n=n+size(useCells,1);
                 % end
                  
end
  % else
      %  sprintf('couldnt find file %d %s',fnum,afile{fnum})
   % end
end


save('JLH_Development','exptdata','unitdata','-v7.3');