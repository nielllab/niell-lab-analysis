function [] = NRKO_lev3_ppc_convol_group_peristim

global info
global outputDir
output  = 'lev3_tuningstats_group_peristim';
filename = fullfile(outputDir,output,output);
mkdir(fullfile(outputDir,output))
load(filename);
%
rate_stat.unitinfo(rate_stat.unitinfo(:,4)==2,4) = 3;

%% plot the different layer 4 cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
gt = {'N2BKO','wt','N2AKO'};

figure
cnt = 0;
for iStim = 1:2
    for k = 0:1
        for iLayer = 3:5
            cors = {'r', 'k', 'g'};            
            df = [];
            cnt=cnt+1;
            subplot(4,3,cnt)
            mn = []; sm = [];
            for iType = [1 0 2]
                hold on
                if k==0
                    sl =  rate_stat.unitinfo(:,4)==iLayer & rate_stat.unitinfo(:,5)==k & rate_stat.unitinfo(:,1)==iType & rate_stat.unitinfo(:,6)==iStim;
                elseif k==1
                    sl =  rate_stat.unitinfo(:,4)==iLayer & rate_stat.unitinfo(:,5)==k & rate_stat.unitinfo(:,1)==iType & rate_stat.unitinfo(:,6)==iStim;
                end                    
                if sum(sl)==0, continue,end
                dat = rate_stat.osi(sl);
                dat(dat==0) = NaN;
                dat(~isfinite(dat)) = NaN;

                [A,B,C] = unique(rate_stat.animalinfo(sl));
                cfg.semiweighted_adaptive = 'yes';
                stt = Lur_statistics(cfg, dat, C, rate_stat.unitinfo(sl,1));
                mn = stt.means_semiweighted;
                sm = stt.sems_semiweighted;         
                h = errorbar(iType,mn,sm,[cors{iType+1} 'o'])
                df(iType+1) = sum(sl);
                                n(iType+1) = stt.nIndividuals;

            end
            yax = get(gca, 'YLim');
            yax(1) = 0;
            set(gca,'YLim', yax);
            xlim([-0.5 2.5])
            nt = {'pyr', 'int'};
            stims = {'drift', 'bar'};
            layers = {'nan', 'nan','superficial', 'granular', 'deep'};
            H = title(sprintf('%s, %s, %s, %d(%d) %d(%d) %d(%d)',  layers{iLayer}, stims{iStim}, nt{k+1},df(1),n(1), df(2),n(2), df(3), n(3)));
            set(H,'FontSize',6)
        end
    end
end
%%
figure
cnt = 0;
for iStim = 1:2
    for k = 0:1
        for iLayer = 3:5
            cors = {'r', 'k', 'g'};            
            df = [];
            cnt=cnt+1;
            subplot(4,3,cnt)
            mn = []; sm = [];
            for iType = [1 0 2]
                hold on
                if k==0
                    sl =  rate_stat.unitinfo(:,4)==iLayer & rate_stat.unitinfo(:,2)==k & rate_stat.unitinfo(:,1)==iType & rate_stat.unitinfo(:,6)==iStim;
                elseif k==1
                    sl =  rate_stat.unitinfo(:,4)==iLayer & rate_stat.unitinfo(:,2)==k & rate_stat.unitinfo(:,1)==iType & rate_stat.unitinfo(:,6)==iStim;
                end                    
                if sum(sl)==0, continue,end
                dat = rate_stat.osi(sl);
                dat(dat==0) = NaN;
                dat(~isfinite(dat)) = NaN;

                [A,B,C] = unique(rate_stat.animalinfo(sl));
                cfg.semiweighted_adaptive = 'yes';
                stt = Lur_statistics(cfg, dat, C, rate_stat.unitinfo(sl,1));
                mn = stt.means_semiweighted;
                sm = stt.sems_semiweighted;         
                h = errorbar(iType,mn,sm,[cors{iType+1} 'o'])
                df(iType+1) = sum(sl);
                                n(iType+1) = stt.nIndividuals;

            end
            yax = get(gca, 'YLim');
            yax(1) = 0;
            set(gca,'YLim', yax);
            xlim([-0.5 2.5])
            set(gca, 'XTickLabels', gt)
            nt = {'non-NR5A', 'NR5A'};
            stims = {'drift', 'bar'};
            layers = {'nan', 'nan','superficial', 'granular', 'deep'};
            H = title(sprintf('%s, %s, %s, %d(%d) %d(%d) %d(%d)',  layers{iLayer}, stims{iStim}, nt{k+1}, df(1),n(1), df(2),n(2), df(3), n(3)));
            set(H,'FontSize',6)
        end
    end
end
%%

gt = {'N2BKO','wt','N2AKO'};

figure
cnt = 0;
for iStim = 1:2
    for k = 0:1
        for iLayer = 3:5
            cors = {'r', 'k', 'g'};            
            df = [];
            cnt=cnt+1;
            subplot(4,3,cnt)
            mn = []; sm = [];
            for iType = [1 0 2]
                hold on
                if k==0
                    sl =  rate_stat.unitinfo(:,4)==iLayer & rate_stat.unitinfo(:,5)==k & rate_stat.unitinfo(:,1)==iType & rate_stat.unitinfo(:,6)==iStim;
                elseif k==1
                    sl =  rate_stat.unitinfo(:,4)==iLayer & rate_stat.unitinfo(:,5)==k & rate_stat.unitinfo(:,1)==iType & rate_stat.unitinfo(:,6)==iStim;
                end                    
                if sum(sl)==0, continue,end
                dat = rate_stat.stim_mod(sl);
                dat(dat==0) = NaN;
                dat(~isfinite(dat)) = NaN;

                [A,B,C] = unique(rate_stat.animalinfo(sl));
                cfg.semiweighted_adaptive = 'yes';
                stt = Lur_statistics(cfg, dat, C, rate_stat.unitinfo(sl,1));
                mn = stt.means_semiweighted;
                sm = stt.sems_semiweighted;         
                h = errorbar(iType,mn,sm,[cors{iType+1} 'o'])
                df(iType+1) = sum(sl);
                                n(iType+1) = stt.nIndividuals;

            end
            ylabel('stim mod')
            yax = get(gca, 'YLim');
            yax(1) = 0;
            set(gca,'YLim', yax);
            xlim([-0.5 2.5])
            nt = {'pyr', 'int'};
            stims = {'drift', 'bar'};
            layers = {'nan', 'nan','superficial', 'granular', 'deep'};
            H = title(sprintf('%s, %s, %s, %d(%d) %d(%d) %d(%d)',  layers{iLayer}, stims{iStim}, nt{k+1},df(1),n(1), df(2),n(2), df(3), n(3)));
            set(H,'FontSize',6)
        end
    end
end
%%
figure
cnt = 0;
for iStim = 1:2
    for k = 0:1
        for iLayer = 3:5
            cors = {'r', 'k', 'g'};            
            df = [];
            cnt=cnt+1;
            subplot(4,3,cnt)
            mn = []; sm = [];
            for iType = [1 0 2]
                hold on
                if k==0
                    sl =  rate_stat.unitinfo(:,4)==iLayer & rate_stat.unitinfo(:,2)==k & rate_stat.unitinfo(:,1)==iType & rate_stat.unitinfo(:,6)==iStim;
                elseif k==1
                    sl =  rate_stat.unitinfo(:,4)==iLayer & rate_stat.unitinfo(:,2)==k & rate_stat.unitinfo(:,1)==iType & rate_stat.unitinfo(:,6)==iStim;
                end                    
                if sum(sl)==0, continue,end
                dat = rate_stat.stim_mod(sl);
                dat(dat==0) = NaN;
                dat(~isfinite(dat)) = NaN;

                [A,B,C] = unique(rate_stat.animalinfo(sl));
                cfg.semiweighted_adaptive = 'yes';
                stt = Lur_statistics(cfg, dat, C, rate_stat.unitinfo(sl,1));
                mn = stt.means_semiweighted;
                sm = stt.sems_semiweighted;         
                h = errorbar(iType,mn,sm,[cors{iType+1} 'o'])
                df(iType+1) = sum(sl);
                                n(iType+1) = stt.nIndividuals;

            end
            yax = get(gca, 'YLim');
            yax(1) = 0;
            ylabel('stim mod')
            set(gca,'YLim', yax);
            xlim([-0.5 2.5])
            set(gca, 'XTickLabels', gt)
            nt = {'non-NR5A', 'NR5A'};
            stims = {'drift', 'bar'};
            layers = {'nan', 'nan','superficial', 'granular', 'deep'};
            H = title(sprintf('%s, %s, %s, %d(%d) %d(%d) %d(%d)',  layers{iLayer}, stims{iStim}, nt{k+1}, df(1),n(1), df(2),n(2), df(3), n(3)));
            set(H,'FontSize',6)
        end
    end
end
%%

gt = {'N2BKO','wt','N2AKO'};

figure
cnt = 0;
for iStim = 1:2
    for k = 0:1
        for iLayer = 3:5
            cors = {'r', 'k', 'g'};            
            df = [];
            cnt=cnt+1;
            subplot(4,3,cnt)
            mn = []; sm = [];
            for iType = [1 0 2]
                hold on
                if k==0
                    sl =  rate_stat.unitinfo(:,4)==iLayer & rate_stat.unitinfo(:,5)==k & rate_stat.unitinfo(:,1)==iType & rate_stat.unitinfo(:,6)==iStim;
                elseif k==1
                    sl =  rate_stat.unitinfo(:,4)==iLayer & rate_stat.unitinfo(:,5)==k & rate_stat.unitinfo(:,1)==iType & rate_stat.unitinfo(:,6)==iStim;
                end                    
                if sum(sl)==0, continue,end
                dat = rate_stat.rate(sl,1,2);
                dat(dat==0) = NaN;
                dat(~isfinite(dat)) = NaN;

                [A,B,C] = unique(rate_stat.animalinfo(sl));
                cfg.semiweighted_adaptive = 'yes';
                stt = Lur_statistics(cfg, dat, C, rate_stat.unitinfo(sl,1));
                mn = stt.means_semiweighted;
                sm = stt.sems_semiweighted;         
                h = errorbar(iType,mn,sm,[cors{iType+1} 'o'])
                df(iType+1) = sum(sl);
                                n(iType+1) = stt.nIndividuals;

            end
            ylabel('stim mod')
            yax = get(gca, 'YLim');
            yax(1) = 0;
            set(gca,'YLim', yax);
            xlim([-0.5 2.5])
            nt = {'pyr', 'int'};
            stims = {'drift', 'bar'};
            layers = {'nan', 'nan','superficial', 'granular', 'deep'};
            H = title(sprintf('%s, %s, %s, %d(%d) %d(%d) %d(%d)',  layers{iLayer}, stims{iStim}, nt{k+1},df(1),n(1), df(2),n(2), df(3), n(3)));
            set(H,'FontSize',6)
        end
    end
end
%%
figure
cnt = 0;
for iStim = 1:2
    for k = 0:1
        for iLayer = 3:5
            cors = {'r', 'k', 'g'};            
            df = [];
            cnt=cnt+1;
            subplot(4,3,cnt)
            mn = []; sm = [];
            for iType = [1 0 2]
                hold on
                if k==0
                    sl =  rate_stat.unitinfo(:,4)==iLayer & rate_stat.unitinfo(:,2)==k & rate_stat.unitinfo(:,5)==0 & rate_stat.unitinfo(:,1)==iType & rate_stat.unitinfo(:,6)==iStim;
                elseif k==1
                    sl =  rate_stat.unitinfo(:,4)==iLayer & rate_stat.unitinfo(:,2)==k & rate_stat.unitinfo(:,5)==0 & rate_stat.unitinfo(:,1)==iType & rate_stat.unitinfo(:,6)==iStim;
                end                    
                if sum(sl)==0, continue,end
                dat = rate_stat.rate(sl,1,1);
                dat(dat==0) = NaN;
                dat(~isfinite(dat)) = NaN;

                [A,B,C] = unique(rate_stat.animalinfo(sl));
                cfg.semiweighted_adaptive = 'yes';
                stt = Lur_statistics(cfg, dat, C, rate_stat.unitinfo(sl,1));
                mn = stt.means_semiweighted;
                sm = stt.sems_semiweighted;         
                h = errorbar(iType,mn,sm,[cors{iType+1} 'o'])
                df(iType+1) = sum(sl);
                                n(iType+1) = stt.nIndividuals;

            end
            yax = get(gca, 'YLim');
            yax(1) = 0;
            ylabel('stim mod')
        %    set(gca,'YLim', yax);
            xlim([-0.5 2.5])
            set(gca, 'XTickLabels', gt)
            nt = {'non-NR5A', 'NR5A'};
            stims = {'drift', 'bar'};
            layers = {'nan', 'nan','superficial', 'granular', 'deep'};
            H = title(sprintf('%s, %s, %s, %d(%d) %d(%d) %d(%d)',  layers{iLayer}, stims{iStim}, nt{k+1}, df(1),n(1), df(2),n(2), df(3), n(3)));
            set(H,'FontSize',6)
        end
    end
end
%%
gt = {'N2BKO','wt','N2AKO'};

figure
for iLatency = 1:2
    figure
    cnt = 0;
    for iStim = 1:2
        for k = 0:1
            for iLayer = 3:5
                cors = {'r', 'k', 'g'};            
                df = [];
                cnt=cnt+1;
                subplot(4,3,cnt)
                mn = []; sm = [];
                for iType = [1 0 2]
                    hold on
                    if k==0
                        sl =  rate_stat.unitinfo(:,4)==iLayer & rate_stat.unitinfo(:,5)==k & rate_stat.unitinfo(:,1)==iType & rate_stat.unitinfo(:,6)==iStim;
                    elseif k==1
                        sl =  rate_stat.unitinfo(:,4)==iLayer & rate_stat.unitinfo(:,5)==k & rate_stat.unitinfo(:,1)==iType & rate_stat.unitinfo(:,6)==iStim;
                    end                    
                    if sum(sl)==0, continue,end

                    dat = rate_stat.loco_mod(sl,iLatency);
                    dat(dat==0) = NaN;
                    dat(~isfinite(dat)) = NaN;
                    [A,B,C] = unique(rate_stat.animalinfo(sl));
                    cfg.semiweighted_adaptive = 'yes';
                    stt = Lur_statistics(cfg, dat, C, rate_stat.unitinfo(sl,1));
                    mn = stt.means_semiweighted;
                    sm = stt.sems_semiweighted;         
                    h = errorbar(iType,mn,sm,[cors{iType+1} 'o'])
                    df(iType+1) = sum(sl);
                    n(iType+1) = stt.nIndividuals;
                end
                ylabel('locomod')
             %   yax = get(gca, 'YLim');
             %   yax(1) = 0;
             %   set(gca,'YLim', yax);
                xlim([-0.5 2.5])
                nt = {'pyr', 'int'};
                latencies = {'prestim', 'stim'};
                stims = {'drift', 'bar'};
                layers = {'nan', 'nan','superficial', 'granular', 'deep'};
                H = title(sprintf('%s %s, %s, %s, %d(%d) %d(%d) %d(%d)',  latencies{iLatency},layers{iLayer}, stims{iStim}, nt{k+1}, df(1),n(1), df(2),n(2), df(3), n(3)));
                set(H,'FontSize',6)
            end
        end
    end
end
%%
for iLatency = 1:2
    figure
    cnt = 0;
    for iStim = 1:2
        for k = 0:1
            for iLayer = 3:5
                try
                    cors = {'r', 'k', 'g'};            
                    df = [];
                    cnt=cnt+1;
                    subplot(4,3,cnt)
                    mn = []; sm = [];
                    for iType = [1 0 2]
                        hold on
                        if k==0
                            sl =  rate_stat.unitinfo(:,4)==iLayer & rate_stat.unitinfo(:,2)==k & rate_stat.unitinfo(:,1)==iType & rate_stat.unitinfo(:,6)==iStim;
                        elseif k==1
                            sl =  rate_stat.unitinfo(:,4)==iLayer & rate_stat.unitinfo(:,2)==k & rate_stat.unitinfo(:,1)==iType & rate_stat.unitinfo(:,6)==iStim;
                        end                    
                        if sum(sl)==0, continue,end
                        dat = rate_stat.loco_mod(sl,iLatency);
                        dat(dat==0) = NaN;
                        dat(~isfinite(dat)) = NaN;

                        [A,B,C] = unique(rate_stat.animalinfo(sl));
                        cfg.semiweighted_adaptive = 'yes';
                        stt = Lur_statistics(cfg, dat, C, rate_stat.unitinfo(sl,1));
                        mn = stt.means_semiweighted;
                        sm = stt.sems_semiweighted;         
                        h = errorbar(iType,mn,sm,[cors{iType+1} 'o'])
                        df(iType+1) = sum(sl);
                                        n(iType+1) = stt.nIndividuals;

                    end
                    yax = get(gca, 'YLim');
                    ylabel('mod by locomotion')
                    %yax(1) = 0;
                    set(gca,'YLim', yax);
                    xlim([-0.5 2.5])
                    set(gca, 'XTickLabels', gt)
                    nt = {'non-NR5A', 'NR5A'};
                    stims = {'drift', 'bar'};
                    layers = {'nan', 'nan','superficial', 'granular', 'deep'};
                    H = title(sprintf('%s %s, %s, %s, %d(%d) %d(%d) %d(%d)',  latencies{iLatency}, layers{iLayer}, stims{iStim}, nt{k+1}, df(1),n(1), df(2),n(2), df(3), n(3)));
                    set(H,'FontSize',6)                    
                end
            end
        end
    end
end