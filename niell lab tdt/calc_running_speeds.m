    
%[fname pname] =uiputfile('*.ps'); psfilename=fullfile(pname,fname);  %%% get ps filename
%psfilename = 'c:/test.ps';   %%% default location
%if exist(psfilename,'file')==2;delete(psfilename);end %%% check for previous file
clear all

 [fname pname] =uiputfile('*.ps'); psfilename=fullfile(pname,fname);  %%% get ps filename
%psfilename = 'c:/test.ps';   %%% default location
if exist(psfilename,'file')==2;delete(psfilename);end

pname = uigetdir('C:\data\tdt tanks','block data')
delims = strfind(pname,'\');
selected_path = pname(1 :delims(length(delims))-1)
Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
Block_Name = pname(delims(length(delims))+1 :length(pname))

nChan = input('number of chans (only need 1): ');
movement = 1;

tic

flags =struct('lfpTseries',1,'lfpSpectra',1,'mouseOn',1);

[afname, apname] = uigetfile('*.mat','analysis data');
afile = fullfile(apname,afname);
load(afile);
use_afile=1;

tdtData= getTDTdata(Tank_Name, Block_Name, 1:nChan, flags);


       
    if movement
        hold on
        tsamp = tdtData.mouseT;
        vsmooth = tdtData.mouseV;
        
       figure
        plot(tsamp,vsmooth,'g');
               
        set(gcf, 'PaperPositionMode', 'auto');
        print('-dpsc',psfilename,'-append');
    end
   
    
   
    %%%%

    
    v_interp = interp1(tsamp,vsmooth);
    
figure
hist(vsmooth,1:2:30)

figure
hist(vsmooth,1:1:30)
title '1cm/s bins'
    %%%%
    
   locomotion.speeds = v_interp;
   locomotion.mouseT = tsamp;
   locomotion.mouseV = vsmooth;
   
    



save(afile,'locomotion','-append');
 

 
 