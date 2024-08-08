%% Figure 5
load('SegEx')
filterWindow = [5 40];
filterOrder = 4;
Fs = 1000;
evaluation_angle = pi;
win_angle = 0;
tol = 0.2;
times = 2000+fixOn+200:2050; %spontaneous window minus 200 ms from fixaton onset
[rows,cols,~] = size( thisLFP ); channels = rows*cols;
[X,Y] = meshgrid( 1:cols, 1:rows );

% Filter and calculate phase
filtLFP = bandpass_filter(thisLFP,filterWindow(1),filterWindow(2),filterOrder,Fs);
figure
plot(squeeze(nanmean(nanmean(filtLFP(:,:,times)))),'k','linewidth',1.5)
axis off
hold on
plot([20 70],[-65 -65],'k','linewidth',1.5)
plot([20 20],[-65 -40],'k','linewidth',1.5)



spikeTimes = [];
for j = 1:length(mySpikes)
    theseSpikes = [];
    for k = 1:size(mySpikes{j},1)
        if mySpikes{j}{k,2} < 5
            theseSpikes = [theseSpikes mySpikes{j}{k,1}' + 2000]; %Put spikes into absolute time instead of relative to stim onset
        end
    end
    if length(find(theseSpikes > times(1) & theseSpikes < times(end))) > 0
        spikeTimes(j,:) = histc(theseSpikes(theseSpikes > times(1) & theseSpikes < times(end))-times(1),1:length(times));
    else
        spikeTimes(j,:) = zeros(1,length(times));
    end
end
chanList = find(sum(spikeTimes,2) > 4);
%     chanList = [2 3 5 6 7 9];
plotRaster(spikeTimes(chanList,:),1:length(times),[0 0 0])
axis off


%
ShuffData1 = [];
ShuffData2 = [];
DSIData = [];
load('TwixMotifTuning.mat')
fractionT = length(find(motifLength > zPrctile))./length(motifLength);
DSIData = allZ(2,:);
ShuffData1 = allZ(1,:);
ShuffData2 = shuffZ(2,:);
load('WolfieMotifTuning.mat')
fractionW = length(find(motifLength > zPrctile))./length(motifLength);
DSIData = [DSIData allZ(2,:)];
ShuffData1 = [ShuffData1 allZ(1,:)];
ShuffData2 = [ShuffData2 shuffZ(2,:)];
[h1 p1] = signrank(DSIData,ShuffData1)
[h2 p2] = signrank(DSIData,ShuffData2)
bins = linspace(-2,10,15);

count = 1;
for ii = 1:2
    if ii == 1
        load('TwixMotifTuning')
        motifEx = [14];
    else
        load('WolfieMotifTuning')
        motifEx = [26];
    end
    for i = 1:length(tuningBin)
        ChanBins(i) = length(tuningBin{i});
    end
    for kk = 1:length(motifEx)
        jj = motifEx(kk);
        if motifLength(jj) > zPrctile(jj)
            figure
            [p] = polar(linspace(0,2*pi,length(motifRate(jj,:))+1),[motifRate(jj,:) motifRate(jj,1)]);
            hold on
            p.Color = 'k';
            p.LineWidth = 1.5;
            [p1] = polar(linspace(0,2*pi,length(motifRate(jj,:))+1),[motifRate(jj,:)-((motifSTD(jj,:))) motifRate(jj,1)-(motifSTD(jj,1))]);
            [p2] = polar(linspace(0,2*pi,length(motifRate(jj,:))+1),[motifRate(jj,:)+((motifSTD(jj,:))) motifRate(jj,1)+(motifSTD(jj,1))]);
            p1.Color = 'k';
            p2.Color = 'k';
            p1.LineStyle = ':';
            p2.LineStyle = ':';
            [r] = polar(linspace(0,2*pi,length(permMean(jj,:))+1),[permMean(jj,:) permMean(jj,1)]);
            r.LineWidth = 1;
            r.Color = 'r';
            [r1] = polar(linspace(0,2*pi,length(permMean(jj,:))+1),[permMean(jj,:)-permCI(jj,:) permMean(jj,1)-permCI(jj,1)]);
            [r2] = polar(linspace(0,2*pi,length(permMean(jj,:))+1),[permMean(jj,:)+permCI(jj,:) permMean(jj,1)+permCI(jj,1)]);
            r1.Color = 'r';
            r2.Color = 'r';
            r1.LineStyle = ':';
            r2.LineStyle = ':';
            rho = motifRate(jj,:);
            theta = linspace(0,2*pi,length(motifRate(jj,:))+1);
            [x y] = pol2cart(theta(1:length(theta)-1),rho);
            thisAngle = cart2pol(sum(x),sum(y));
            [q1] = quiver(0,0,cos(thisAngle),sin(thisAngle),motifLength(jj));
            q1.LineWidth = 1.5;
            q1.Color = 'k';
            rho = permMean(jj,:);
            theta = linspace(0,2*pi,length(permMean(jj,:))+1);
            [x y] = pol2cart(theta(1:length(theta)-1),rho);
            thisAngle = cart2pol(sum(x),sum(y));
            thisLength = sqrt((sum(x).^2)+(sum(y).^2));
            [q2] = quiver(0,0,cos(thisAngle),sin(thisAngle),thisLength);
            q2.LineWidth = 1.5;
            q2.Color = 'r';
            xlim([-2 2])
            ylim([-2 2])
            title(['Title ' num2str(jj) 'Motif Z = ' num2str(allZ(2,jj))])
            print(['PanelB' num2str(count) '.pdf'],'-painters','-dpdf')
            count = count+1;
        end
    end
end


figure
histogram(DSIData,bins,'facecolor','b','normalization','probability','linewidth',1.5)
hold on
histogram(ShuffData2,bins,'facecolor','r','normalization','probability','linewidth',1.5)
set(gca,'fontsize',16,'linewidth',1.5)
xlabel('Tuning Z-Score')
ylabel('Fraction')
box off
title('Tuning Shuffle')
legend('Data','Shuffle')

