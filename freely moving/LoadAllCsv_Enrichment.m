savePDF=1; dbstop if error

fileList=[] ;fileListR=[] ;fileListL=[] ; TSfileList=[]; %finds all files w/top.csv in the name
for i=1:length(pname)
    fileList = [fileList; dir([pname{i} '*resnet50_Top*.csv'])];  
    
end

 name =fileList.name;

%%

for j=1:length(fileList)
    if savePDF
        psfilename = 'C:\analysisPSControl.ps';
        if exist(psfilename,'file')==2;delete(psfilename);end
    end
    clear path 
       
    fname=fullfile(fileList(j).folder,fileList(j).name);

    sname = split(fname,'_');
    ani = sname{2}(end-4:end);
    sessionnum = sname{4}(end);
    date = sname{5};
    clipnum = sname{6}(~isletter(sname{6}));
    
%     sname = split(fname,'_');
%     ani = sname{1}(end-4:end);
%     sessionnum = sname{3}(end);
%     date = sname{4};
%     clipnum = sname{5}(~isletter(sname{5}));
       
    Data(j).ani= {ani};
    Data(j).sessionnum = {sessionnum}; Data(j).date = {date};
    Data(j).clipnum = {clipnum};
    Data(j).Data = csvread(fname,3,0);
   
    aligned = alignHeadEnrichment(fname,4,0,psfilename,.99, .95)
    Data(j).mouse_xy=aligned.mouse_xy;
    Data(j).mouseV=((aligned.mouseSp)*27)./30;
    Data(j).theta=aligned.theta;
    Data(j).dTheta=aligned.dTheta;
    Data(j).cricketxy=aligned.crick_xy;
    Data(j).cricketV=((aligned.crickSp)*27)./30;
    Data(j).range=((aligned.range)*27)./30;
    Data(j).az=aligned.az;
    Data(j).cricketP=aligned.crick_p;
    
   
    if savePDF
        pSname='T:\PreyCaptureAnalysis\Data\ControlAnalysis\';
        C={ani, date, sessionnum, clipnum};
        filen=sprintf('%s%s%s%s',ani,date,sessionnum,clipnum,'.pdf')
        pdfilename=fullfile(pSname,filen)
        dos(['ps2pdf ' psfilename ' ' pdfilename]);
        delete(psfilename);
    end 
end
afilename=sprintf('%s',ani,'_control','.mat')
save(fullfile(pSname, afilename))

