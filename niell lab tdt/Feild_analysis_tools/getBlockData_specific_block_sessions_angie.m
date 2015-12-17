

clear all

dbstop if error
n=0;

%D:\Angie_analysis\DOI_experiments\09_23_15\analysis_092315.mat 
analysisPath = 'D:\Angie_analysis\DOI_experiments\'

for dataset =1:1
   
    if dataset==1
    afile = {'09_23_15\analysis_092315.mat',...
   '11_19_15\analysis_111915.mat',...
   '12_02_15\analysis_120215.mat'}
   
    elseif dataset==2
        
    analysisPath = 'D:\Jen_analysis\NR5A_Pinping\Jen_NR5A_analysis_files\analysis_files\N2B\'

    afile = {};

    else dataset==3
        analysisPath = 'D:\Jen_analysis\NR5A_Pinping\Jen_NR5A_analysis_files\analysis_files\N2A\'

    afile = {}
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
         
%         for blocknum = 3:4
            Block_Name
            B = strncmpi(Block_Name,'bar',1) ;
            
                 if sum(B)>1 ;
            blocknum = listdlg('ListString',Block_Name,'SelectionMode','single');
                 else
            blocknum=find(B==1);
                 end
                  
            if sum(B)~=0
                      
            Block_Name{blocknum};
            
            flags = struct('visStim',1,'lfpTseries',1,'lfpSpectra',1,'mouseOn',1);
            % flags = struct('visStim',1);
            
            ch = unique(useCells(:,1))';
          
            ch_lfp=1:64;
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
                    %unitdata{cell_n+n}.responsive = resp(cell_n);
                    unitdata{cell_n+n}.layer=layer(cell_n);
                    %unitdata{cell_n+n}.inhibitory = inh(cell_n);
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


save('Angie_DOI_bar_LFP_unit','exptdata','unitdata','-v7.3');