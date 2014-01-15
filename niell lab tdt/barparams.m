        function barparams(parameter,plist,ncol)
        if ~exist('ncol','var') | isempty(ncol)
            ncol=5;
        end
        meds = [];
        errs = [];
        for i = 1:ncol;
            meds=[meds nanmean(parameter(plist{i}))];
            errs = [errs stderr(parameter(plist{i}))];
        end
        %barweb(meds,errs);
        bar(meds);
        hold on
        errorbar(meds,errs,'k.');
        xlim([0 ncol+1])
        end
