[fps pps] = uigetfile('*.ps', 'ps file');
[fpdf ppdf] = uiputfile('*.pdf','pdf file');


    dos(['ps2pdf "' fullfile(pps,fps) '" "' fullfile(ppdf,fpdf)  '"'] )

    
    