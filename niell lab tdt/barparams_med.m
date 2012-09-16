        function barparams_med(parameter,plist)
        
        meds = [];
        errs = [];
        for i = 1:5;
            meds=[meds nanmedian(parameter(plist{i}))];
            errs = [errs stderr(parameter(plist{i}))];
        end
        bar(meds);
        hold on
        errorbar(meds,errs,'k.');
        xlim([0 6])
        end
