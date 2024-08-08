%% Create a pinwheel
theta = linspace(0,2*pi,5000);
theta(end) = [];
shift = pi/2;
x = linspace(0+shift,8*pi+shift,length(theta));
y = (sin(x).*(1-(1/sqrt(2))))+1;
newTheta = randsample(theta,length(theta),true,y);
newTheta = sort(newTheta);

rho = linspace(1,99,99);
side = 60;
screen = zeros(200,200);
[X Y] = meshgrid(1:size(screen,2),1:size(screen,1));
x0 = size(screen,2)/2;
y0 = size(screen,1)/2;

for i = 1:length(theta)
    for j = 1:length(rho)
        [x y] = pol2cart(theta(i),rho(j));
        screen(y0-round(y),x0-round(x)) = newTheta(i);
    end
end
screen = screen(y0-floor(side/2)+1:y0+floor(side/2),x0-floor(side/2)+1:x0+floor(side/2));
figure
imagesc(screen)
axis square
map = colorcet( 'C2' ); colormap( circshift( map, [ 28, 0 ] ) )
axis off
colorbar

bins = linspace(0,2*pi,101);
bins(end) = [];
binStep = bins(2)-bins(1);
newScreen = screen;
for i = 1:length(bins)
    newScreen(screen >= bins(i) & screen < bins(i)+binStep) = bins(i);
end
oldScreen = newScreen;
meanBin = length(newScreen(:))./length(bins);
binCount = [];
for i = 1:length(bins)
    binCount(i) = length(find(newScreen(:) == bins(i)));
end
shuff = [1 2];
linScreen = newScreen(:);
binOff = binCount-meanBin;
figure
plot(binOff)
while max(abs(binOff)) > 1
    for i = 1:length(binCount)        
        if binOff(i) > 0
            theseBins = find(linScreen == bins(i));
            randBins = theseBins(randperm(length(theseBins),abs(binOff(i))));
            for j = 1:length(randBins)
                binShift = circshift(bins,shuff(randi(2,1)));
                linScreen(randBins(j)) = binShift(i);
            end
        elseif binOff(i) < 0
            bins1 = circshift(bins,-shuff(randi(2,1)));
            bins2 = circshift(bins,shuff(randi(2,1)));
            theseBins = find(linScreen == bins1(i) | linScreen == bins2(i));
            randBins = theseBins(randperm(length(theseBins),abs(binOff(i))));
            linScreen(randBins) = bins(i);                
        end
        binCount = [];
        for k = 1:length(bins)
            binCount(k) = length(find(linScreen == bins(k)));
        end
    end
    binOff = binCount-meanBin;
    plot(binOff)
    pause(0.1)
end
close
low = find(binOff < 0);
high = find(binOff > 0);
for i = 1:length(low)
    theseBins = find(linScreen == bins(low(i)));
    thoseBins = find(linScreen == bins(high(i)));
    linScreen(thoseBins(1)) = bins(low(i));
end
newScreen = reshape(linScreen,60,60);

figure
imagesc(newScreen)
axis square
map = colorcet( 'C2' ); colormap( circshift( map, [ 28, 0 ] ) )
axis off
colorbar

figure
histogram(newScreen(:),bins)

%% Flip pinwheel to create full 4mm map
test = newScreen;
test2 = rot90(test,2);
test2 = fliplr(test2);
test3 = rot90(test2,2);
test4 = rot90(test2,2);
test4 = flipud(test4);

allTest = zeros(size(test,1)*2,size(test,2)*2);
allTest(1:size(test,1),1:size(test,2)) = test;
allTest(size(test,1)+1:size(test,1)*2,1:size(test,2)) = test2;
allTest(1:size(test,1),size(test,2)+1:size(test,2)*2) = test3;
allTest(size(test,1)+1:size(test,1)*2,size(test,2)+1:size(test,2)*2) = test4;

test = allTest;
test2 = rot90(test,2);
test2 = fliplr(test2);
test3 = rot90(test2,2);
test4 = rot90(test2,2);
test4 = flipud(test4);

allTest = zeros(size(test,1)*2,size(test,2)*2);
allTest(1:size(test,1),1:size(test,2)) = test;
allTest(size(test,1)+1:size(test,1)*2,1:size(test,2)) = test2;
allTest(1:size(test,1),size(test,2)+1:size(test,2)*2) = test3;
allTest(size(test,1)+1:size(test,1)*2,size(test,2)+1:size(test,2)*2) = test4;

test = allTest;
test2 = rot90(test,2);
test2 = fliplr(test2);
test3 = rot90(test2,2);
test4 = rot90(test2,2);
test4 = flipud(test4);

allTest = zeros(size(test,1)*2,size(test,2)*2);
allTest(1:size(test,1),1:size(test,2)) = test;
allTest(size(test,1)+1:size(test,1)*2,1:size(test,2)) = test2;
allTest(1:size(test,1),size(test,2)+1:size(test,2)*2) = test3;
allTest(size(test,1)+1:size(test,1)*2,size(test,2)+1:size(test,2)*2) = test4;

test = allTest;
test2 = rot90(test,2);
test2 = fliplr(test2);
test3 = rot90(test2,2);
test4 = rot90(test2,2);
test4 = flipud(test4);

allTest = zeros(size(test,1)*2,size(test,2)*2);
allTest(1:size(test,1),1:size(test,2)) = test;
allTest(size(test,1)+1:size(test,1)*2,1:size(test,2)) = test2;
allTest(1:size(test,1),size(test,2)+1:size(test,2)*2) = test3;
allTest(size(test,1)+1:size(test,1)*2,size(test,2)+1:size(test,2)*2) = test4;
screen = allTest(1:600,1:600);
allTest = screen;
figure
imagesc(screen)
map = colorcet( 'C2' ); colormap( circshift( map, [ 28, 0 ] ) )
axis square
c = colorbar;
c.Label.String = 'Preferred Direction (rad)';
c.TickLabels = {'0','\pi','2 \pi'};
c.Ticks = [0 pi 2*pi];
set(gca,'fontsize',16,'linewidth',1.5)
set(gca,'xtick',150:150:600,'xticklabel',1:1:4)
set(gca,'ytick',1:150:451,'yticklabel',4:-1:1)
xlabel('Distance (mm)')
ylabel('Distance (mm)')
