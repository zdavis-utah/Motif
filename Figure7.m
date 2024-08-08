%% Figure 7
load('EvokedResults')

figure
[b] = boxchart([zeros(1,length(find(sigLength < 1))) zeros(1,length(find(sigLength < 1)))+1 zeros(1,length(find(sigLength > 1)))+2],[abs(randAngle(randperm(length(randAngle),length(find(sigLength < 1))))) abs(sigAngle(sigLength < 1)) abs(sigAngle(sigLength > 1))],'boxfacecolor','k','notch','on','linewidth',1.5,'markercolor','k','markerstyle','.');
set(gca,'xtick',[0 1 2],'xticklabel',{'Random Motifs','Weak Modulation','Strong Modulation'},'fontsize',16,'linewidth',1.5,'ytick',[0 pi/2 pi],'yticklabel',{'0','\pi/2','\pi'})
ylabel('Perf-Motif Angle Difference (rad)')
ylim([0 pi])
[h p] = ranksum(abs(randAngle(randperm(length(randAngle),length(find(sigLength < 1))))),abs(sigAngle(sigLength > 1)))
[h p] = ranksum(abs(sigAngle(sigLength < 1)),abs(sigAngle(sigLength > 1)))


motifTrack = [motifTrack motifTrack];
figure
scatter([zeros(1,length(sigUpMod))+1 zeros(1,length(sigDownMod))+2],[sigUpMod sigDownMod],20,'filled','markerfacecolor',[.4 .4 .4])
hold on
scatter([1 2],[nanmean(sigUpMod) nanmean(sigDownMod)],100,'filled','markerfacecolor','k')
for i = 1:length(sigUpMod)
    plot([1 2],[sigUpMod(i) sigDownMod(i)],'color',[0.4 0.4 0.4],'linewidth',motifTrack(i)/75)
end
plot([1 2],[nanmean(sigUpMod) nanmean(sigDownMod)],'k','linewidth',2)
set(gca,'fontsize',15,'linewidth',1.5,'xtick',[1 2],'xticklabel',{'Motif Aligned','Motif Opposed'},'ytick',[0.5 1 1.5],'yticklabel',[-50 0 +50])
xlim([0.5 2.5])
ylim([0.25 1.75])
ylabel('Perf. Mod. (%)')
motifTrack = motifTrack./max(motifTrack);
[nanmean(sigUpMod.*motifTrack) nanmean(sigDownMod.*motifTrack)]
[nanstd(sigUpMod.*motifTrack)./sqrt(length(sigUpMod)-1) nanstd(sigDownMod.*motifTrack)./sqrt(length(sigDownMod)-1)]
[h p] = signrank(sigUpMod.*motifTrack,sigDownMod.*motifTrack)

figure
hold on
bins = linspace(-0.75,0.75,17);
index1 = (sigUpMod-sigDownMod)./(sigUpMod+sigDownMod);
index2 = (sigRandUp-sigRandDown)./(sigRandUp+sigRandDown);
histogram(index1,bins,'normalization','probability','linewidth',1.5)
histogram(index2,bins,'normalization','probability','linewidth',1.5)
set(gca,'fontsize',16,'linewidth',1.5,'xtick',-1:0.5:1,'xticklabel',-100:50:100)
xlabel('Performance Modulation (%)')
ylabel('Fraction of Motifs')
legend('Motifs','Shuffle')



