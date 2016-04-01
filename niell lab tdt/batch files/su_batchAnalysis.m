%clear all

% BatchEphys

for f = 1:length(files)
    
    f
    
    analysisFile = [files(f).path files(f).analysisfile];
    load(analysisFile);
   
    clusterFile = [files(f).path files(f).clusterfile];
    load(clusterFile);
    
    pdfname = [analysisFile(1:end-4) '.pdf'];
    
    % white noise analysis
    for i = 1:length(files(f).blockWn)
        [y blocknum] = ismember(files(f).blockWn{i},Block_Name);
        
        if y==0
        error ('block name does not exist, check that bloack name entered correctly')
        end
        
        movieFile = 'C:\Users\lab\Desktop\movieFiles\cortex\wn_cortex_012alpha1_5hzLg30Hz.mat';
        load(movieFile);
        
        if ismember (Block_Name(blocknum),files(f).postBlocks);
            post=1;
        else
            post=0;
        end
        
        %movietype is entered after blocknum input and equals 1 for contrast
        %modualted noise movie
        noise_analysis(clusterFile,analysisFile,pdfname,movieFile,Block_Name{blocknum},blocknum,1,1,i);
        close all
    end
    % drift analysis
   
    for i=1:length(files(f).blockDrift);
        [y blocknum] = ismember(files(f).blockDrift{i},Block_Name);
        
        if y==0
        error ('block name does not exist, check that bloack name entered correctly')
        end
        
        if ismember (Block_Name(blocknum),files(f).postBlocks);
            post=1;
        else
            post=0;
        end
        
        drift_analysis_move(clusterFile,analysisFile,pdfname,Block_Name{blocknum},blocknum,post)
        close all
    end
    % Optogenetic tagging 
    
    for i=1:length(files(f).blockPinp);
        [y blocknum] = ismember(files(f).blockPinp{i},Block_Name);
        
        if y==0
        error ('block name does not exist, check that bloack name entered correctly')
        end
        
        
        nchan=files(f).nchan;
       
        OptoTag(clusterFile, analysisFile, pdfname,Block_Name{blocknum},blocknum,nchan,i);
        close all
    end
    %% adding layer info to analysis file
%     
     appending_drift_layer(analysisFile, files(f).tip_loc_1,files(f).tip_loc_2,files(f).angle);
%    

end


%  condition(f) = files(f).condition;
%     sex(f) = files(f).sex;
%     age(f) = files(f).age;