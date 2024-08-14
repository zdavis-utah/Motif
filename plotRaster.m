function [] = plotRaster(spikeMat, tVec,color)
fig = figure;
fig.Color = [1 1 1];
%set(gcf,'position',[100 500 1000 100])
%axis off
hold on;
for trialCount = 1:size(spikeMat,1)
    %chStatus = strcat('Working on channel #',num2str(trialCount))
    spikeTimes = find(spikeMat(trialCount,:) == 1);
    spikePos = tVec(spikeTimes);
    for spikeCount = 1:length(spikePos)
        plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [trialCount-0.4 trialCount+0.4],'color',color,'linewidth',2);
    end
end
ylim([0 size(spikeMat, 1)+1]);