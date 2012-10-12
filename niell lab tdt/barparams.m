        function barparams(parameter,plist)
        
        meds = [];
        errs = [];
        for i = 1:5;
            meds=[meds nanmean(parameter(plist{i}))];
            errs = [errs stderr(parameter(plist{i}))];
        end
        %barweb(meds,errs);
        bar(meds);
        hold on
        errorbar(meds,errs,'k.');
        xlim([0 6])
        end
