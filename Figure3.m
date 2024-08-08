%% Figure 3

load( "F3_data.mat")
kvals = {'Random','BPN'};
% F3b -- similarity distributions
figure; hold on
for k = 1:length(kvals) 
    plot( sim_dist(k).x, sim_dist(k).f , 'linewidth',2)
end
legend( kvals)
xlabel( 'Similarity') 

% F3c -- ratio plot
figure; hold on
for k = 1:length(kvals) 
    if ~isempty( dist_ratio(k).ratio)
    plot( k, dist_ratio(k).ratio , 'k.', 'markersize', 50)
    plot( [k-0.2,k+0.2], [mean(dist_ratio(k).ratio), mean(dist_ratio(k).ratio)], '-k', 'linewidth',2)
    else 
        plot(k, 0, 'k.')
    end
end
set(gca, 'XTick', 1:length(kvals), 'XTickLabels', kvals)
xlim( [0, length(kvals) + 1])
ylabel( 'Cluster Distance Ratio'); xlabel( 'Clustering K-Value') 


% Panel D
load('ModelMotifData125_Good.mat')
theta = tuneBins;
motifAngle = []; motifLength = [];
ii = 3;
[x y] = pol2cart(theta(1:end-1),(meanRate(ii,1:end-1)./tuneRate(1:end-1))-0.8);
[motifAngle(ii) motifLength(ii)] = cart2pol(sum(x),sum(y));

figure
[p1] = polarplot(tuneBins, [meanRate(ii,1:end-1)./tuneRate(1:end-1) meanRate(ii,1)./tuneRate(end)]);
p1.LineWidth = 1.5;
p1.Color = 'k';
hold on
[p2] = polarplot(tuneBins,[meanRate(ii,1:end-1)./tuneRate(1:end-1) meanRate(ii,1)./tuneRate(end)]-meanSTD(ii,:)./tuneRate);
p2.LineWidth = 1.5;
p2.Color = 'k';
p2.LineStyle = ':';
[p3] = polarplot(tuneBins,[meanRate(ii,1:end-1)./tuneRate(1:end-1) meanRate(ii,1)./tuneRate(end)]+meanSTD(ii,:)./tuneRate);
p3.LineWidth = 1.5;
p3.Color = 'k';
p3.LineStyle = ':';
[p4] = polarplot(tuneBins,[meanRand(ii,1:end-1) meanRand(ii,1)]);
p4.LineWidth = 1.5;
p4.Color = 'r';
[p5] = polarplot(tuneBins,[meanRand(ii,1:end-1) meanRand(ii,1)]-randSTD(ii,:));
p5.LineWidth = 1.5;
p5.Color = 'r';
p5.LineStyle = ':';
[p6] = polarplot(tuneBins,[meanRand(ii,1:end-1) meanRand(ii,1)]+randSTD(ii,:));
p6.LineWidth = 1.5;
p6.Color = 'r';
p6.LineStyle = ':';
polarplot([motifAngle(ii) motifAngle(ii)],[0.8 motifLength(ii)+0.8],'k','linewidth',2)
rlim([0.8 1.1])


load('ModelMotifData0_Good.mat')
theta = tuneBins;
motifAngle = []; motifLength = [];
ii = 2;
[x y] = pol2cart(theta(1:end-1),(meanRate(ii,1:end-1)./tuneRate(1:end-1))-0.8);
[motifAngle(ii) motifLength(ii)] = cart2pol(sum(x),sum(y));


figure
[p1] = polarplot(tuneBins, [meanRate(ii,1:end-1)./tuneRate(1:end-1) meanRate(ii,1)./tuneRate(end)]);
p1.LineWidth = 1.5;
p1.Color = 'k';
hold on
[p2] = polarplot(tuneBins,[meanRate(ii,1:end-1)./tuneRate(1:end-1) meanRate(ii,1)./tuneRate(end)]-meanSTD(ii,:)./tuneRate);
p2.LineWidth = 1.5;
p2.Color = 'k';
p2.LineStyle = ':';
[p3] = polarplot(tuneBins,[meanRate(ii,1:end-1)./tuneRate(1:end-1) meanRate(ii,1)./tuneRate(end)]+meanSTD(ii,:)./tuneRate);
p3.LineWidth = 1.5;
p3.Color = 'k';
p3.LineStyle = ':';
[p4] = polarplot(tuneBins,[meanRand(ii,1:end-1) meanRand(ii,1)]);
p4.LineWidth = 1.5;
p4.Color = 'r';
[p5] = polarplot(tuneBins,[meanRand(ii,1:end-1) meanRand(ii,1)]-randSTD(ii,:));
p5.LineWidth = 1.5;
p5.Color = 'r';
p5.LineStyle = ':';
[p6] = polarplot(tuneBins,[meanRand(ii,1:end-1) meanRand(ii,1)]+randSTD(ii,:));
p6.LineWidth = 1.5;
p6.Color = 'r';
p6.LineStyle = ':';
polarplot([motifAngle(ii) motifAngle(ii)],[0.8 motifLength(ii)+0.8],'k','linewidth',2)
rlim([0.8 1.1])


% Panel E
load('HistVars2')
figure
histogram(zLength,linspace(-3,8,20),'normalization','probability','linewidth',1.5)
hold on
histogram(zRand,linspace(-3,8,20),'normalization','probability','linewidth',1.5)
set(gca,'fontsize',16,'linewidth',1.5)
box off
xlabel('Tuning Z-Score')
ylabel('Fraction of Time')
legend('Data','Shuffle')

