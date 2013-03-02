function [] = GroupPsychoPhysPlot(groupPsy)

%cd '/Users/alangordon/mounts/w5/alan/perceptMnemonic/behav_data'
printit = 0;
plotMean = 1;
plotIndSubjects = 0;

plotScatter = 0;
plotFitLine = 0;
plotFitCurve = 1;

plotCurve = 1;
groupMinN = 0;

maxGraphArarySize = 20;

fn = fieldnames(groupPsy(1).dat);

NSubs = size(groupPsy,2);

if plotMean
    m = ceil(sqrt(maxGraphArarySize+1));
    NPlots = NSubs+1;
else
    m = ceil(sqrt(maxGraphArarySize));
    NPlots = NSubs;
end

subplotDims = [5 4];

for f =1:length(fn)
    
    
        page = 0;
        for i = 1:length(groupPsy)
            
            if plotIndSubjects
                i2 = 1 + mod(i-1,maxGraphArarySize);
            else
                i2 = i;
            end
                
            if (i2==1);
                
                if (printit && (i~=1))
                    h = gcf;
                    print(h, '-depsc', [fn{f} 'part' num2str(page)])
                end
                page = page+1;
                figure('Color', 'w', 'Position', [100 100 750 750], 'PaperPositionMode', 'auto');
                iGoodSubs = 0;
                hold on;
            end
            
            if ~isempty(groupPsy(i).dat)
                thisPlotDat = groupPsy(i).dat.(fn{f});
                
                plotColor='k';
                
                if plotScatter
                    if plotIndSubjects
                        subplot(subplotDims(1),subplotDims(2),i2)
                        scatter(thisPlotDat.x, thisPlotDat.y)
                        
                        if plotFitLine
                            xCoords = [min(thisPlotDat.x), max(thisPlotDat.x)];
                            yCoords = [thisPlotDat.betas(1) + thisPlotDat.betas(2) * min(thisPlotDat.x) ...
                                thisPlotDat.betas(1) + thisPlotDat.betas(2) * max(thisPlotDat.x)];
                            line(xCoords, yCoords);
                            
                            subMeans.betas(i,:) = thisPlotDat.betas;
                        end
                        
                        if plotFitCurve
                            hold on;
                            plot(thisPlotDat.x, thisPlotDat.fitCurve, '-r', 'LineWidth', 2);
                        end
                    end
                end
                
                
                if plotCurve
                    %if groupPsy(i).goodSub
                    plotColor='b';
                    iGoodSubs = iGoodSubs+1;
                    subMeans.(fn{f})(iGoodSubs,:) = thisPlotDat.mean;
                    %end
                    
                    if plotIndSubjects
                         subplot(subplotDims(1),subplotDims(2),i2)
                        errorbar(thisPlotDat.xVals, thisPlotDat.mean, thisPlotDat.SE,  thisPlotDat.marker, 'Color', plotColor, 'MarkerFaceColor', plotColor, 'MarkerSize', 4);
                        
                        if isfield(thisPlotDat, 'g_fc')
                            hold on
                            plot(thisPlotDat.g_coh_signed, thisPlotDat.g_fc, plotColor);
                            subMeans.betas(iGoodSubs,:) = thisPlotDat.g_fc;
                        end
                    end
                end
                
                fullXLabel = arrayfun(@(x) {num2str(x)}, thisPlotDat.ticks.x);
                
                set(gca, 'XLim', [min(thisPlotDat.ticks.x)  max(thisPlotDat.ticks.x)], 'XTick', thisPlotDat.ticks.x, 'XTickLabel', fullXLabel, ...
                    'YLim', thisPlotDat.ticks.y, 'TickDir', 'out', 'Color', 'w','box','off');
                
                grid on;
                xlabel(thisPlotDat.xlabel);
                ylabel(thisPlotDat.ylabel);
                
            end
        end
    
        
        %%
        if plotScatter
            if plotMean
                
                i = i+1;
                if plotIndSubjects
                    subplot(m,m,i)
                end
                
                meanBetas = mean(subMeans.betas);
                
                xCoords = [min(thisPlotDat.ticks.x), max(thisPlotDat.ticks.x)];
                yCoords = [meanBetas(1) + meanBetas(2) * min(thisPlotDat.ticks.x) ...
                    meanBetas(1) + meanBetas(2) * max(thisPlotDat.ticks.x)];
                
                line(xCoords, yCoords);
                
                set(gca, 'XLim', [min(thisPlotDat.ticks.x) - .5  max(thisPlotDat.ticks.x) + .5], 'XTick', thisPlotDat.ticks.x, 'XTickLabel', makeTickLabel(thisPlotDat.ticks.x,1), ...
                    'YLim', thisPlotDat.ticks.y, 'TickDir', 'out', 'Color','w','box','off');
                grid on;
                xlabel(thisPlotDat.xlabel);
                ylabel(thisPlotDat.ylabel);
            end
        end
        
        %%
        if plotCurve
            if plotMean
                groupMean = nanmean(subMeans.(fn{f}));
                groupSE = nanstd(subMeans.(fn{f})) ./ sqrt(sum(~isnan(subMeans.(fn{f}))));
                
                groupN = sum(~isnan(subMeans.(fn{f})));
                removeLowNs_h = double(groupN>groupMinN);
                removeLowNs = removeLowNs_h;
                removeLowNs(removeLowNs_h==0) = NaN;
                
                groupMean = groupMean .* removeLowNs;
                
                groupSE = groupSE .* removeLowNs;
                %groupSE = thisPlotDat.groupAnovaSE .* removeLowNs;
                
                i2 = i2+1;
                if plotIndSubjects
                subplot(subplotDims(1),subplotDims(2),i2)                    
                end
                
                errorbar(thisPlotDat.xVals, groupMean, groupSE, '-o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
                
                if isfield(thisPlotDat, 'g_fc')
                    hold on
                    
                    g_fc = mean(subMeans.betas);
                    plot(thisPlotDat.g_coh_signed, g_fc, 'r');
                end
                
                fullXLabel = arrayfun(@(x) {num2str(x)}, thisPlotDat.ticks.x);
                
                set(gca, 'XLim', [min(thisPlotDat.ticks.x) - .5  max(thisPlotDat.ticks.x) + .5], 'XTick', thisPlotDat.ticks.x, 'XTickLabel', fullXLabel, ...
                    'YLim', thisPlotDat.ticks.y, 'TickDir', 'out', 'Color','w','box','off');
                grid off;
                xlabel(thisPlotDat.xlabel);
                ylabel(thisPlotDat.ylabel);
            end
        end
        
        if printit
            h = gcf;
            print(h, '-depsc', [fn{f} 'part' num2str(page)])
        end
    end
end
