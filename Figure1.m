%% Figure 1
% Make Connection Probability Map Figure
% Dependencies: circ-dist from P. Berens, CircStat: A Matlab Toolbox for 
% Circular Statistics, Journal of Statistical Software, Volume 31, Issue 
% 10, 2009. colorcet.m from Peter Kovesi. Good Colour Maps: How to Design 
% Them. Model connectivity generating code MakeModelConnections.m and
% MakeTuningMap.m. Simulations performed using NETSIM: https://github.com/mullerlab/NETSIM

close all; clear all;
Xe = [300:-1:0 1:300];
Ye = [300:-1:0 1:300];
numCon = 1000;
for type = 1:2
    if type == 2
        load('allTest')
        a = 0.1; % Amplitude of Gaussian
        k = 1.25; % Concentration of von Miesis.
        b = 1;
        c = 72; % Gaussian Sigma: 72 = 600 µm
        fracLoc = 0.9; %Local Weight
        fracClust = 0.1; %TuningWeight
        saveID = {'E','F'};
        load('ConT')
    else
        load('allTest')
        a = 0.1; % Amplitude of Gaussian
        k = 0;
        b = 1;
        c = 72;
        fracLoc = 0.9;
        fracClust = 0.1;
        saveID = {'B','C'};
        load('ConNT')
    end
    
    thisLocE = [300 300];
    thatLocE = [262 262]; %Location of example neuron to form connections
    thisRho = allTest(thisLocE(1),thisLocE(2));
    dist = [];
    spatProbE = [];
    tuningProbE = [];
    for j = 1:size(allTest,1)
        for i = 1:size(allTest,2)
            dist = sqrt((Xe(j)-Xe(thisLocE(1))).^2 + (Ye(i)-Ye(thisLocE(2))).^2);
            spatProbE(i,j) = a*exp(-((dist).^2)/(2*c.^2)); % Gaussian
            tuningDist = circ_dist(thisRho,allTest(i,j));
            tuningProbE(i,j) = b*(exp(k*cos(tuningDist))./(2*pi)); % Von Miesis
        end
    end
    remapLoc = [thatLocE(1) thatLocE(2)];
    if remapLoc(1) <= thisLocE(1)
        shiftX1 = Xe(remapLoc(1))+1:size(spatProbE,1);
        shiftX2 = 1:(size(spatProbE,1)-length(shiftX1));
    else
        shiftX1 = size(spatProbE,1)-Xe(remapLoc(1)):size(spatProbE,1);
        shiftX2 = 1:(size(spatProbE,1)-length(shiftX1));
    end
    if remapLoc(2) <= thisLocE(1)
        shiftY1 = Ye(remapLoc(2))+1:size(spatProbE,1);
        shiftY2 = 1:(size(spatProbE,2)-length(shiftY1));
    else
        shiftY1 = size(spatProbE,2)-Ye(remapLoc(2)):size(spatProbE,1);
        shiftY2 = 1:(size(spatProbE,2)-length(shiftY1));
    end
    thisSpatProb = spatProbE([shiftY1 shiftY2],[shiftX1 shiftX2]);    
    tuningProbE = zeros(size(allTest,1),size(allTest,2));
    for jj = 1:size(allTest,1)
        for ii = 1:size(allTest,2)
            tuningDist = circ_dist(allTest(thatLocE(1),thatLocE(2)),allTest(ii,jj));
            tuningProbE(ii,jj) = b*(exp(k*cos(tuningDist))./(2*pi)); % Von Miesis
        end
    end
    tuningProbE = tuningProbE./max(max(tuningProbE));
    tuningProbE = tuningProbE*.1;
    thisSpatProb = (thisSpatProb*fracLoc)+(fracClust*tuningProbE);
     
    figure
    imagesc(thisSpatProb)
    axis square
    c = colorbar;
    colormap('hot')
    c.Label.String = 'Connection Probability';
    xlabel('Distance (mm)')
    ylabel('Distance (mm)')
    c.Ticks = 0:.05:.1;
    caxis([0 .1])
    set(gca,'xtick',[0 120 240 360 480 600],'xticklabel',0:5,'ytick',[1 120 240 360 480 600],'yticklabel',5:-1:0,'fontsize',16,'linewidth',1.5)
       
    bins = unique(allTest);
    if type == 1
       for i = 1:600
           for j = 1:600
               allTest(i,j) = bins(randi(length(bins),1));
           end
       end
    end
    
    figure
    imagesc(allTest)
    map = colorcet( 'C2' ); map = circshift(map, [28,0]); colormap( map);
    axis square
    c = colorbar;
    c.Label.String = 'Unit Tuning (rad)';
    c.Ticks = 0:pi/2:2*pi;
    c.TickLabels = {'0','\pi/2','\pi','3\pi/2','2\pi'};
    xlabel('Distance (mm)')
    ylabel('Distance (mm)')
    set(gca,'xtick',[0 120 240 360 480 600],'xticklabel',0:5,'ytick',[1 120 240 360 480 600],'yticklabel',5:-1:0,'fontsize',16,'linewidth',1.5) 
    hold on
    [r c] = ind2sub(600,thisCon(1:800));
    scatter(c,r,10,'filled','markerfacecolor','k')
    xlim([1 600])
    ylim([1 600])
    set(gca,'xtick',[0 120 240 360 480 600],'xticklabel',0:5,'ytick',[-600:120:0],'yticklabel',0:5,'fontsize',16)
    axis square
    axis off
        xlabel('Distance (mm)')
    ylabel('Distance (mm)')
end

% Panel G
load('ConnEx2')
figure
imagesc(allTest)
map = colorcet( 'C2' ); map = circshift(map, [28,0]); colormap( map);
axis image
axis off
xlim([1 60])
ylim([1 60])

% Panel H
figure
hold on
linPhase = allTest(:);
phaseBins = linspace(0,2*pi,11);
colorSteps = linspace(1,size(map,1),11);
colorSteps = round(colorSteps);
unitPhase = [];
for ii = 1:1000
    unitPhase(ii) = linPhase(unitIdx(ii));
    phaseBin = histc(unitPhase(ii),phaseBins);
    thisPhase = find(phaseBin == 1);
    [r c] = ind2sub(600,connEx(ii,1:800));
    scatter(r,c,100,'filled','markerfacecolor',map(colorSteps(thisPhase),:))
end
axis image
axis off
xlim([1 60])
ylim([600-60 600])
alpha(0.5)

% PanelI
linTest = allTest(:);
bins = linspace(0,pi,11);
conHist = [];
for i = 1:length(randIdx)
    [c r] = ind2sub(600,randIdx(i));
    thisTune = allTest(r,c);
    thoseCons = linTest(randUnits(i,1:800));
    conHist(i,:) = histc(abs(circ_dist(thisTune,thoseCons)),bins);
end
conHist = conHist./sum(nanmean(conHist));
figure
errorbar(bins(1:end-1),nanmean(conHist(:,1:end-1)),nanstd(conHist(:,1:end-1)),'linewidth',1.5,'color','k')
xlabel('Tuning Distance (rad)')
xlim([-0.1 pi])
ylim([0.05 0.16])
ylabel('Connection Probability')
set(gca,'fontsize',16,'linewidth',1.5,'xtick',bins([1 6 11]),'xticklabel',{'0','\pi/2','\pi'},'ytick',[0.05 0.1 0.15])
box off


