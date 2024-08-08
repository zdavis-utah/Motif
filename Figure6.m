%%Figure 7

% Panel A
load('EvokedResults')
timeBins = -200:10:500;
figure
hold on
[b1 b2] = boundedline(1:size(unitHistAlignedGood,2),nanmean(unitHistAlignedGood),nanstd(unitHistAlignedGood)./sqrt(size(unitHistAlignedGood,1)-1),'b');
[r1 r2] = boundedline(1:size(unitHistOpposedGood,2),nanmean(unitHistOpposedGood),nanstd(unitHistOpposedGood)./sqrt(size(unitHistOpposedGood,1)-1),'r');
b1.LineWidth = 1.5;
r1.LineWidth = 1.5;
alpha(0.5);
xlim([0 70])
set(gca,'xtick',0:10:70,'xticklabel',timeBins(1:10:71),'fontsize',16,'linewidth',1.5)
xlabel('Time from target onset (ms)')
ylabel('Normalized firing rate')
legend([b1 r1],{'Motif Aligned Targets','Motif Opposed Targets'})
% print -painters -dpdf PanelA.pdf
[h1 p1] = signrank(nanmean(unitHistAlignedGood(:,28:43)),nanmean(unitHistOpposedGood(:,28:43)))
[nanmean(nanmean(unitHistAlignedGood(:,28:43))) nanmean(nanmean(unitHistOpposedGood(:,28:43)))]
[nanstd(nanmean(unitHistAlignedGood(:,28:43)))./sqrt(size(unitHistAlignedGood,1)-1) nanstd(nanmean(unitHistOpposedGood(:,28:43)))./sqrt(size(unitHistOpposedGood,1)-1)]

% Panel B
figure
scatter([zeros(1,length(unitAlignedGood))+1 zeros(1,length(unitOpposedGood))+2],[unitAlignedGood unitOpposedGood],30,'filled','markerfacecolor','k')
hold on
scatter([1 2],[nanmean(unitAlignedGood) nanmean(unitOpposedGood)],120,'filled','markerfacecolor','k')
for i = 1:length(unitAlignedGood)
    plot([1 2],[unitAlignedGood(i) unitOpposedGood(i)],'k','linewidth',0.5)
end
plot([1 2],[nanmean(unitAlignedGood) nanmean(unitOpposedGood)],'k','linewidth',3)
[nanstd(unitAlignedGood)./sqrt(size(unitAlignedGood,2)-1) nanstd(unitOpposedGood)./sqrt(size(unitOpposedGood,2)-1)];

set(gca,'fontsize',15,'linewidth',1.5,'xtick',[1 2],'xticklabel',{'Motif Aligned','Motif Opposed'})
xlim([0.5 2.5])
ylim([0 2])
ylabel('Evoked gain modulation (norm.)')
[h p] = signrank(unitAlignedGood,unitOpposedGood)

% Panel C
index1 = (unitAlignedGood-unitOpposedGood)./(unitAlignedGood+unitOpposedGood);
index2 = (unitAlignedGoodRand-unitOpposedGoodRand)./(unitAlignedGoodRand+unitOpposedGoodRand);
figure
histogram(index1,linspace(-1,1,15),'linewidth',1.5,'normalization','probability')
hold on
histogram(index2,linspace(-1,1,15),'linewidth',1.5,'normalization','probability')
xlim([-1 1])
set(gca,'fontsize',16,'linewidth',1.5,'xtick',[-1 -0.5 0 0.5 1])
box off
xlabel('Modulation Index')
ylabel('Fraction of units')
legend('Data','Shuffle')
