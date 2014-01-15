        function barparams_med(parameter,plist,ncol)
                if ~exist('ncol','var') | isempty(ncol)
            ncol=5;
                end
        
        meds = [];
        errs = [];
        for i = 1:ncol;
            meds=[meds nanmedian(parameter(plist{i}))];
            errs = [errs stderr(parameter(plist{i}))];
        end
        bar(meds);
        hold on
        errorbar(meds,errs,'k.');
        xlim([0 ncol+1])
        end
