cd('\\nadata\snlrdata\Zac\Figures\Manuscript\Motif\Model')
clear all;
load('allTest6')
tuneBins = sort(unique(allTest(:)));
kScan = [2];
for scan = 1:length(kScan)
    Xe = [300:-1:0 1:300];
    Ye = [300:-1:0 1:300];
    numCon = 1000;
    a = 0.1; % Amplitude of Gaussian
    k = kScan(scan); % Concentration of von Miesis. Vary linspace(0,4,5)?
    b = 1;
    c = 72; % Gaussian Sigma: 75 = 500 µm so 80% = 400 µm
    tuneWidth = 60;
    fracLoc = .9; %vary ratio [1 .9 .8 .5 0]
    fracClust = .1;
    thisLocE = [300 300]; %Location of neuron to form connections
    thisRho = allTest(thisLocE(1),thisLocE(2));
    dist = [];
    spatProbE = [];
    tuningProbE = [];
    for i = 1:size(allTest,1)
        for j = 1:size(allTest,2)
            dist = sqrt((Xe(j)-Xe(thisLocE(1))).^2 + (Ye(i)-Ye(thisLocE(2))).^2);
            spatProbE(i,j) = a*exp(-((dist).^2)/(2*c.^2)); % Gaussian
        end
    end
    
    % Inhibitory connection probability
    Xi = [150:-1:0 1:150];
    Yi = [150:-1:0 1:150];
    IMap = allTest;
    IMap(2:2:600,:) = [];
    IMap(:,2:2:600) = [];
    spatProbI = spatProbE;
    spatProbI(2:2:600,:) = [];
    spatProbI(:,2:2:600) = [];
    thisLocI = ceil(thisLocE./2);
    count = 1;
    count2 = 360001;
    unitCon = zeros(450000,numCon);
    conLocE = [];
    conWeight = [];
    conCount = 1;
    for i = 1:600
        for j = 1:600
            remapLoc = [i j];
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
            tuningProbE = zeros(tuneWidth,tuneWidth);
            for ii = 1:tuneWidth
                for jj = 1:tuneWidth
                    tuningDist = circ_dist(allTest(i,j),allTest(ii,jj));
                    tuningProbE(ii,jj) = b*(exp(k*cos(tuningDist))./(2*pi)); % Von Miesis
                end
            end
            test = tuningProbE; test2 = rot90(test,2); test2 = fliplr(test2); test3 = rot90(test2,2); test4 = rot90(test2,2); test4 = flipud(test4);
            
            fullTuningProb = zeros(size(test,1)*2,size(test,2)*2);
            fullTuningProb(1:size(test,1),1:size(test,2)) = test;
            fullTuningProb(size(test,1)+1:size(test,1)*2,1:size(test,2)) = test2;
            fullTuningProb(1:size(test,1),size(test,2)+1:size(test,2)*2) = test3;
            fullTuningProb(size(test,1)+1:size(test,1)*2,size(test,2)+1:size(test,2)*2) = test4;
            
            test = fullTuningProb; test2 = rot90(test,2); test2 = fliplr(test2); test3 = rot90(test2,2); test4 = rot90(test2,2); test4 = flipud(test4);
            
            fullTuningProb = zeros(size(test,1)*2,size(test,2)*2);
            fullTuningProb(1:size(test,1),1:size(test,2)) = test;
            fullTuningProb(size(test,1)+1:size(test,1)*2,1:size(test,2)) = test2;
            fullTuningProb(1:size(test,1),size(test,2)+1:size(test,2)*2) = test3;
            fullTuningProb(size(test,1)+1:size(test,1)*2,size(test,2)+1:size(test,2)*2) = test4;
            
            test = fullTuningProb; test2 = rot90(test,2); test2 = fliplr(test2); test3 = rot90(test2,2); test4 = rot90(test2,2); test4 = flipud(test4);
            
            fullTuningProb = zeros(size(test,1)*2,size(test,2)*2);
            fullTuningProb(1:size(test,1),1:size(test,2)) = test;
            fullTuningProb(size(test,1)+1:size(test,1)*2,1:size(test,2)) = test2;
            fullTuningProb(1:size(test,1),size(test,2)+1:size(test,2)*2) = test3;
            fullTuningProb(size(test,1)+1:size(test,1)*2,size(test,2)+1:size(test,2)*2) = test4;
            
            test = fullTuningProb; test2 = rot90(test,2); test2 = fliplr(test2); test3 = rot90(test2,2); test4 = rot90(test2,2); test4 = flipud(test4);
            
            fullTuningProb = zeros(size(test,1)*2,size(test,2)*2);
            fullTuningProb(1:size(test,1),1:size(test,2)) = test;
            fullTuningProb(size(test,1)+1:size(test,1)*2,1:size(test,2)) = test2;
            fullTuningProb(1:size(test,1),size(test,2)+1:size(test,2)*2) = test3;
            fullTuningProb(size(test,1)+1:size(test,1)*2,size(test,2)+1:size(test,2)*2) = test4;
            tuningProbE = fullTuningProb(1:600,1:600);
            tuningProbE = tuningProbE./max(max(tuningProbE));
            tuningProbE = tuningProbE*0.1;
            

            justSpatProb = thisSpatProb;
            thisSpatProb = (thisSpatProb*fracLoc)+(fracClust*tuningProbE);
            justSpatProb(justSpatProb > 0.05) = thisSpatProb(justSpatProb > 0.05);

            linProb1 = thisSpatProb(:);
            linTune = allTest(:);
            idxE = 1:length(linProb1);
            numConE = numCon*.8; %80 percent excitatory connections
            numConI = numCon*.2;

            conLoc = [];
            while length(conLoc) < numConE
                conLoc = [conLoc randsample(idxE,numConE,true,linProb1)]; % Matlab doesn't do weighted random sampling without replacement, so need to correct double connections
                conLoc = unique(conLoc);
            end
            randRmv = randperm(length(conLoc),length(conLoc)-numConE);
            conLoc(randRmv) = [];

            % Generate Inhibitory Connections
            conLocI = conLoc(randperm(length(conLoc),numConI));
            
            [r c] = ind2sub([600 600],conLocI);
            ri = ceil(r./2);
            ci = ceil(c./2);
            conLocI = sub2ind([300 300],ri,ci)+360000;             
            unitCon(count,:) = [conLoc sort(conLocI)];
            count = count+1;
                       
            if ismember(i,2:2:size(allTest,1)) && ismember(j,2:2:size(allTest,2))
                thisSpatProb = justSpatProb(:);
                linProb = thisSpatProb(:);
                idxE = 1:length(linProb);
                numConE = numCon*.8; %80 percent excitatory connections
                numConI = numCon*.2;
                conLoc = [];
                while length(conLoc) < numConE
                    conLoc = [conLoc randsample(idxE,numConE,true,linProb)]; % Matlab doesn't do weighted random sampling without replacement, so need to correct double connections
                    conLoc = unique(conLoc);
                end
                randRmv = randperm(length(conLoc),length(conLoc)-numConE);
                conLoc(randRmv) = [];
                
                randCon = conLoc(randperm(length(conLoc),numConI));
                [r c] = ind2sub([600 600],randCon);
                ri = ceil(r./2);
                ci = ceil(c./2);
                conLocI = sub2ind([300 300],ri,ci)+360000;
                unitCon(count2,:) = [conLoc sort(conLocI)];
                count2 = count2+1;
            end
            if mod(count,1000) == 0
                display(['Percent Complete = ' num2str(count/360000)])
            end
        end
    end
    unitCon = unitCon';
    fid = fopen(['C' num2str(kScan(scan)) 'S600R91-test6.dat'], "wb");
    fwrite(fid, 450000, "uint32");
    fwrite(fid, numCon, "uint32");
    fwrite(fid, unitCon, "uint32");
    fclose(fid);
end
