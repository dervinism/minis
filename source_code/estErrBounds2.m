function [AmpsMean, RTsMean, TwoDsMean, AmpsMedian, RTsMedian, TwoDsMedian, AmpsMin, RTsMin, TwoDsMin, AmpsMax, RTsMax, TwoDsMax, AmpsPrct, RTsPrct, TwoDsPrct,...
    boundAmps, boundRTs, boundTwoDs, boundAmpsLinear, boundRTsLinear, boundTwoDsLinear, AmpsMeanBottom, AmpsMedianBottom, AmpsMinBottom, AmpsMaxBottom,...
	  AmpsPrctBottom, boundAmpsBottom, boundAmpsLinearBottom, AmpsMeanMid, AmpsMedianMid, AmpsMinMid, AmpsMaxMid, AmpsPrctMid, boundAmpsMid, boundAmpsLinearMid,...
	  AmpsMeanTop, AmpsMedianTop, AmpsMinTop, AmpsMaxTop, AmpsPrctTop, boundAmpsTop, boundAmpsLinearTop, sLength, slFactor, fileNames, fileSweeps, flFactor, pow, powWgh,...
	  estCount, UFest, AmpsLSE, RTsLSE, TwoDsLSE, AmpsLinearLSE, RTsLinearLSE, TwoDsLinearLSE, AmpsLSEBottom, AmpsLinearLSEBottom, AmpsLSEMid, AmpsLinearLSEMid,...
	  AmpsLSETop, AmpsLinearLSETop, SD, optimData, dataF, cntF, F, G, H, A, B, C, FBottom, ABottom, FMid, AMid, FTop, ATop] = estErrBounds2(Amps, AmpsNeg, RTs, RTsNeg,...
    twoDs, twoDsNeg, fileNames, fileSweeps, SD, fileAmps, fileRTs, fileAmpsNeg, fileRTsNeg, sLength, detectionParameters, classificationParameters,...
    sweepcount, expand)



%% Calculate the differences between sweeps and their combinations:
initAmps = zeros(size(Amps));
initAmpsNeg = initAmps;
initRTs = zeros(size(RTs));
initRTsNeg = initRTs;
inittwoDs = struct([]);
inittwoDsNeg = struct([]);
slFactor = max(sLength)./sLength;
for ns = 1:sweepcount
    initAmps(ns,:) = Amps(ns,:)*slFactor(ns);
    initAmpsNeg(ns,:) = AmpsNeg(ns,:)*slFactor(ns);
    initRTs(ns,:) = RTs(ns,:)*slFactor(ns);
    initRTsNeg(ns,:) = RTsNeg(ns,:)*slFactor(ns);
    inittwoDs{ns} = twoDs{ns}*slFactor(ns);
    inittwoDsNeg{ns} = twoDsNeg{ns}*slFactor(ns);
end
optimData.initAmps = initAmps;
optimData.initAmpsNeg = initAmpsNeg;
optimData.initRTs = initRTs;
optimData.initRTsNeg = initRTsNeg;
optimData.initTwoDs = inittwoDs;
optimData.initTwoDsNeg = inittwoDsNeg;
%
flFactor = zeros(1,length(fileSweeps));
for nf = 1:length(fileSweeps)
    if nf == 1
        flFactor(nf) = sum(slFactor(1:fileSweeps(nf)));
    else
        flFactor(nf) = sum(slFactor(fileSweeps(nf-1)+1:fileSweeps(nf-1)+fileSweeps(nf)));
    end
end
flFactor = flFactor/min(flFactor);
[dataF1, dataF2] = plotData(classificationParameters.amplitudeArray, classificationParameters.riseTimeArray, fileAmps, fileRTs,...
    detectionParameters.RTinterval, fileNames, max(flFactor)./flFactor, 'Mini EPSP');
[dataF3, dataF4] = plotData(classificationParameters.amplitudeArray, classificationParameters.riseTimeArray, fileAmpsNeg, fileRTsNeg,...
    detectionParameters.RTinterval, fileNames, max(flFactor)./flFactor, 'Mini IPSP');
dataF = [dataF1 dataF2 dataF3 dataF4];
%
diffAmps = struct([]);
diffAmpsDev = diffAmps;
diffAmpsNeg = diffAmps;
diffAmpsNegDev = diffAmps;
diffRTs = diffAmps;
diffRTsDev = diffAmps;
diffRTsNeg = diffAmps;
diffRTsNegDev = diffAmps;
diffTwoDs = diffAmps;
diffTwoDsDev = diffAmps;
diffTwoDsNeg = diffAmps;
diffTwoDsNegDev = diffAmps;
%
optimAmps = diffAmps;
optimAmpsNeg = diffAmps;
optimRTs = diffAmps;
optimRTsNeg = diffAmps;
optimTwoDs = diffAmps;
optimTwoDsNeg = diffAmps;
%
largestSweep = floor(sweepcount/2);
meanAmpsDiff = zeros(largestSweep,1);
meanAmpsDevDiff = meanAmpsDiff;
meanAmpsNegDiff = meanAmpsDiff;
meanAmpsNegDevDiff = meanAmpsDiff;
meanRTsDiff = meanAmpsDiff;
meanRTsDevDiff = meanAmpsDiff;
meanRTsNegDiff = meanAmpsDiff;
meanRTsNegDevDiff = meanAmpsDiff;
meanTwoDsDiff = meanAmpsDiff;
meanTwoDsDevDiff = meanAmpsDiff;
meanTwoDsNegDiff = meanAmpsDiff;
meanTwoDsNegDevDiff = meanAmpsDiff;
%
medianAmpsDiff = meanAmpsDiff;
medianAmpsDevDiff = meanAmpsDiff;
medianAmpsNegDiff = meanAmpsDiff;
medianAmpsNegDevDiff = meanAmpsDiff;
medianRTsDiff = meanAmpsDiff;
medianRTsDevDiff = meanAmpsDiff;
medianRTsNegDiff = meanAmpsDiff;
medianRTsNegDevDiff = meanAmpsDiff;
medianTwoDsDiff = meanAmpsDiff;
medianTwoDsDevDiff = meanAmpsDiff;
medianTwoDsNegDiff = meanAmpsDiff;
medianTwoDsNegDevDiff = meanAmpsDiff;
%
prctAmpsDiff = meanAmpsDiff;
prctAmpsDevDiff = meanAmpsDiff;
prctAmpsNegDiff = meanAmpsDiff;
prctAmpsNegDevDiff = meanAmpsDiff;
prctRTsDiff = meanAmpsDiff;
prctRTsDevDiff = meanAmpsDiff;
prctRTsNegDiff = meanAmpsDiff;
prctRTsNegDevDiff = meanAmpsDiff;
prctTwoDsDiff = meanAmpsDiff;
prctTwoDsDevDiff = meanAmpsDiff;
prctTwoDsNegDiff = meanAmpsDiff;
prctTwoDsNegDevDiff = meanAmpsDiff;
%
minAmpsDiff = meanAmpsDiff;
minAmpsDevDiff = meanAmpsDiff;
minAmpsNegDiff = meanAmpsDiff;
minAmpsNegDevDiff = meanAmpsDiff;
minRTsDiff = meanAmpsDiff;
minRTsDevDiff = meanAmpsDiff;
minRTsNegDiff = meanAmpsDiff;
minRTsNegDevDiff = meanAmpsDiff;
minTwoDsDiff = meanAmpsDiff;
minTwoDsDevDiff = meanAmpsDiff;
minTwoDsNegDiff = meanAmpsDiff;
minTwoDsNegDevDiff = meanAmpsDiff;
%
maxAmpsDiff = meanAmpsDiff;
maxAmpsDevDiff = meanAmpsDiff;
maxAmpsNegDiff = meanAmpsDiff;
maxAmpsNegDevDiff = meanAmpsDiff;
maxRTsDiff = meanAmpsDiff;
maxRTsDevDiff = meanAmpsDiff;
maxRTsNegDiff = meanAmpsDiff;
maxRTsNegDevDiff = meanAmpsDiff;
maxTwoDsDiff = meanAmpsDiff;
maxTwoDsDevDiff = meanAmpsDiff;
maxTwoDsNegDiff = meanAmpsDiff;
maxTwoDsNegDevDiff = meanAmpsDiff;
%
nCompare = zeros(1,largestSweep);
nSweeps = ceil(sweepcount./(1:largestSweep));
reduceComp = nSweeps - floor(sweepcount./(1:largestSweep));
progress = 100/largestSweep: 100/largestSweep :100;
h = waitbar(0,'Estimating optimisation error bounds - 0%');
waitText1 = 'Estimating optimisation error bounds - ';
waitText2 = '%';
for sweepSize = 1:largestSweep
    nCompare(sweepSize) = nchoosek(nSweeps(sweepSize),2);
    if reduceComp(sweepSize)
        nCompare(sweepSize) = nCompare(sweepSize) - 1;
    end
    Amps = zeros(nSweeps(sweepSize),size(initAmps,2));
    AmpsNeg = Amps;
    RTs = zeros(nSweeps(sweepSize),size(initRTs,2));
    RTsNeg = RTs;
    twoDs = struct([]);
    twoDsNeg = struct([]);
    for newSweep = 1:nSweeps(sweepSize)
        if newSweep ~= nSweeps(sweepSize)
            Amps(newSweep,:) = sum(initAmps((newSweep-1)*sweepSize+1: newSweep*sweepSize,:),1);
            AmpsNeg(newSweep,:) = sum(initAmpsNeg((newSweep-1)*sweepSize+1: newSweep*sweepSize,:),1);
            RTs(newSweep,:) = sum(initRTs((newSweep-1)*sweepSize+1: newSweep*sweepSize,:),1);
            RTsNeg(newSweep,:) = sum(initRTsNeg((newSweep-1)*sweepSize+1: newSweep*sweepSize,:),1);
            twoDs{newSweep} = inittwoDs{(newSweep-1)*sweepSize+1};
            twoDsNeg{newSweep} = inittwoDsNeg{(newSweep-1)*sweepSize+1};
            for iCell = (newSweep-1)*sweepSize+2: newSweep*sweepSize
                twoDs{newSweep} = twoDs{newSweep} + inittwoDs{iCell};
                twoDsNeg{newSweep} = twoDsNeg{newSweep} + inittwoDsNeg{iCell};
            end
        else
            Amps(newSweep,:) = sum(initAmps(end-sweepSize+1: end,:),1);
            AmpsNeg(newSweep,:) = sum(initAmpsNeg(end-sweepSize+1: end,:),1);
            RTs(newSweep,:) = sum(initRTs(end-sweepSize+1: end,:),1);
            RTsNeg(newSweep,:) = sum(initRTsNeg(end-sweepSize+1: end,:),1);
            twoDs{newSweep} = inittwoDs{sweepcount-sweepSize+1};
            twoDsNeg{newSweep} = inittwoDsNeg{sweepcount-sweepSize+1};
            for iCell = sweepcount-sweepSize+2: sweepcount
                twoDs{newSweep} = twoDs{newSweep} + inittwoDs{iCell};
                twoDsNeg{newSweep} = twoDsNeg{newSweep} + inittwoDsNeg{iCell};
            end
        end
    end
    [diffAmps{sweepSize}, diffAmpsDev{sweepSize}, diffAmpsNeg{sweepSize}, diffAmpsNegDev{sweepSize}, diffRTs{sweepSize}, diffRTsDev{sweepSize},...
        diffRTsNeg{sweepSize}, diffRTsNegDev{sweepSize}, diffTwoDs{sweepSize}, diffTwoDsDev{sweepSize}, diffTwoDsNeg{sweepSize},...
        diffTwoDsNegDev{sweepSize}] = compareSweeps(reduceComp(sweepSize), Amps, AmpsNeg, RTs, RTsNeg, twoDs, twoDsNeg);
    
    optimAmps{sweepSize} = Amps;
    optimAmpsNeg{sweepSize} = AmpsNeg;
    optimRTs{sweepSize} = RTs;
    optimRTsNeg{sweepSize} = RTsNeg;
    optimTwoDs{sweepSize} = twoDs;
    optimTwoDsNeg{sweepSize} = twoDsNeg;
    
    meanAmpsDiff(sweepSize) = mean(diffAmps{sweepSize});
    meanAmpsDevDiff(sweepSize) = mean(diffAmpsDev{sweepSize});
    meanAmpsNegDiff(sweepSize) = mean(diffAmpsNeg{sweepSize});
    meanAmpsNegDevDiff(sweepSize) = mean(diffAmpsNegDev{sweepSize});
    meanRTsDiff(sweepSize) = mean(diffRTs{sweepSize});
    meanRTsDevDiff(sweepSize) = mean(diffRTsDev{sweepSize});
    meanRTsNegDiff(sweepSize) = mean(diffRTsNeg{sweepSize});
    meanRTsNegDevDiff(sweepSize) = mean(diffRTsNegDev{sweepSize});
    meanTwoDsDiff(sweepSize) = mean(diffTwoDs{sweepSize});
    meanTwoDsDevDiff(sweepSize) = mean(diffTwoDsDev{sweepSize});
    meanTwoDsNegDiff(sweepSize) = mean(diffTwoDsNeg{sweepSize});
    meanTwoDsNegDevDiff(sweepSize) = mean(diffTwoDsNegDev{sweepSize});
    
    medianAmpsDiff(sweepSize) = median(diffAmps{sweepSize});
    medianAmpsDevDiff(sweepSize) = median(diffAmpsDev{sweepSize});
    medianAmpsNegDiff(sweepSize) = median(diffAmpsNeg{sweepSize});
    medianAmpsNegDevDiff(sweepSize) = median(diffAmpsNegDev{sweepSize});
    medianRTsDiff(sweepSize) = median(diffRTs{sweepSize});
    medianRTsDevDiff(sweepSize) = median(diffRTsDev{sweepSize});
    medianRTsNegDiff(sweepSize) = median(diffRTsNeg{sweepSize});
    medianRTsNegDevDiff(sweepSize) = median(diffRTsNegDev{sweepSize});
    medianTwoDsDiff(sweepSize) = median(diffTwoDs{sweepSize});
    medianTwoDsDevDiff(sweepSize) = median(diffTwoDsDev{sweepSize});
    medianTwoDsNegDiff(sweepSize) = median(diffTwoDsNeg{sweepSize});
    medianTwoDsNegDevDiff(sweepSize) = median(diffTwoDsNegDev{sweepSize});
    
    minAmpsDiff(sweepSize) = min(diffAmps{sweepSize});
    minAmpsDevDiff(sweepSize) = min(diffAmpsDev{sweepSize});
    minAmpsNegDiff(sweepSize) = min(diffAmpsNeg{sweepSize});
    minAmpsNegDevDiff(sweepSize) = min(diffAmpsNegDev{sweepSize});
    minRTsDiff(sweepSize) = min(diffRTs{sweepSize});
    minRTsDevDiff(sweepSize) = min(diffRTsDev{sweepSize});
    minRTsNegDiff(sweepSize) = min(diffRTsNeg{sweepSize});
    minRTsNegDevDiff(sweepSize) = min(diffRTsNegDev{sweepSize});
    minTwoDsDiff(sweepSize) = min(diffTwoDs{sweepSize});
    minTwoDsDevDiff(sweepSize) = min(diffTwoDsDev{sweepSize});
    minTwoDsNegDiff(sweepSize) = min(diffTwoDsNeg{sweepSize});
    minTwoDsNegDevDiff(sweepSize) = min(diffTwoDsNegDev{sweepSize});
    
    maxAmpsDiff(sweepSize) = max(diffAmps{sweepSize});
    maxAmpsDevDiff(sweepSize) = max(diffAmpsDev{sweepSize});
    maxAmpsNegDiff(sweepSize) = max(diffAmpsNeg{sweepSize});
    maxAmpsNegDevDiff(sweepSize) = max(diffAmpsNegDev{sweepSize});
    maxRTsDiff(sweepSize) = max(diffRTs{sweepSize});
    maxRTsDevDiff(sweepSize) = max(diffRTsDev{sweepSize});
    maxRTsNegDiff(sweepSize) = max(diffRTsNeg{sweepSize});
    maxRTsNegDevDiff(sweepSize) = max(diffRTsNegDev{sweepSize});
    maxTwoDsDiff(sweepSize) = max(diffTwoDs{sweepSize});
    maxTwoDsDevDiff(sweepSize) = max(diffTwoDsDev{sweepSize});
    maxTwoDsNegDiff(sweepSize) = max(diffTwoDsNeg{sweepSize});
    maxTwoDsNegDevDiff(sweepSize) = max(diffTwoDsNegDev{sweepSize});
    
    % Estimate the 50th centile of combined EPSP and IPSP data:
    if sweepSize == 1
        Bottom = .5;
        Mid = .9;
        Top = .98;
        if expand
            combAmps = sum(initAmps,1);
            Ampcs = cumsum(combAmps)/sum(combAmps);
            combAmpsNeg = sum(initAmpsNeg,1);
            AmpcsNeg = cumsum(combAmpsNeg)/sum(combAmpsNeg);
            
            initAmpsBottom = initAmps;
            Ampcs = find(Ampcs > Bottom, 1);
            initAmpsBottom(:,1:Ampcs-1) = zeros(size(initAmpsBottom,1), Ampcs-1);
            initAmpsNegBottom = initAmpsNeg;
            AmpcsNeg = find(AmpcsNeg > Bottom, 1);
            initAmpsNegBottom(:,1:AmpcsNeg-1) = zeros(size(initAmpsNegBottom,1), AmpcsNeg-1);
            [diffAmpsBottom, diffAmpsNegBottom] = compareSweepsPartial(reduceComp, initAmpsBottom, initAmpsNegBottom);
            
            initAmpsMid = initAmps;
            Ampcs = find(Ampcs > Mid, 1);
            initAmpsMid(:,1:Ampcs-1) = zeros(size(initAmpsMid,1), Ampcs-1);
            initAmpsNegMid = initAmpsNeg;
            AmpcsNeg = find(AmpcsNeg > Mid, 1);
            initAmpsNegMid(:,1:AmpcsNeg-1) = zeros(size(initAmpsNegMid,1), AmpcsNeg-1);
            [diffAmpsMid, diffAmpsNegMid] = compareSweepsPartial(reduceComp, initAmpsMid, initAmpsNegMid);
            
            initAmpsTop = initAmps;
            Ampcs = find(Ampcs > Top, 1);
            initAmpsTop(:,1:Ampcs-1) = zeros(size(initAmpsTop,1), Ampcs-1);
            initAmpsNegTop = initAmpsNeg;
            AmpcsNeg = find(AmpcsNeg > Top, 1);
            initAmpsNegTop(:,1:AmpcsNeg-1) = zeros(size(initAmpsNegTop,1), AmpcsNeg-1);
            [diffAmpsTop, diffAmpsNegTop] = compareSweepsPartial(reduceComp, initAmpsTop, initAmpsNegTop);
            
            prct = calc50thCentileExpand(diffAmps{1}, diffRTs{1}, diffTwoDs{1}, diffAmpsBottom, diffAmpsMid, diffAmpsTop, 50);
            prctNeg = calc50thCentileExpand(diffAmpsNeg{1}, diffRTsNeg{1}, diffTwoDsNeg{1}, diffAmpsNegBottom, diffAmpsNegMid, diffAmpsNegTop, 50);
        else
            prct = calc50thCentile(diffAmps{1}, diffRTs{1}, diffTwoDs{1}, diffAmpsNegBottom, diffAmpsNegMid, diffAmpsNegTop, 50);
            prctNeg = calc50thCentile(diffAmpsNeg{1}, diffRTsNeg{1}, diffTwoDsNeg{1}, 50);
        end
        optimData.prct = prct;
        optimData.prctNeg = prctNeg;
    end
    prctAmpsDiff(sweepSize) = prctile(diffAmps{sweepSize},prct);
    prctAmpsDevDiff(sweepSize) = prctile(diffAmpsDev{sweepSize},prct);
    prctAmpsNegDiff(sweepSize) = prctile(diffAmpsNeg{sweepSize},prctNeg);
    prctAmpsNegDevDiff(sweepSize) = prctile(diffAmpsNegDev{sweepSize},prctNeg);
    prctRTsDiff(sweepSize) = prctile(diffRTs{sweepSize},prct);
    prctRTsDevDiff(sweepSize) = prctile(diffRTsDev{sweepSize},prct);
    prctRTsNegDiff(sweepSize) = prctile(diffRTsNeg{sweepSize},prctNeg);
    prctRTsNegDevDiff(sweepSize) = prctile(diffRTsNegDev{sweepSize},prctNeg);
    prctTwoDsDiff(sweepSize) = prctile(diffTwoDs{sweepSize},prct);
    prctTwoDsDevDiff(sweepSize) = prctile(diffTwoDsDev{sweepSize},prct);
    prctTwoDsNegDiff(sweepSize) = prctile(diffTwoDsNeg{sweepSize},prctNeg);
    prctTwoDsNegDevDiff(sweepSize) = prctile(diffTwoDsNegDev{sweepSize},prctNeg);
    
    waitText = [waitText1 num2str(progress(sweepSize)) waitText2];
    waitbar(progress(sweepSize)/100, h, waitText);
end
close(h);
X = (1:sweepSize);
AmpsMean = [X; meanAmpsDiff'; meanAmpsDevDiff'; meanAmpsNegDiff'; meanAmpsNegDevDiff'];
RTsMean = [X; meanRTsDiff'; meanRTsDevDiff'; meanRTsNegDiff'; meanRTsNegDevDiff'];
TwoDsMean = [X; meanTwoDsDiff'; meanTwoDsDevDiff'; meanTwoDsNegDiff'; meanTwoDsNegDevDiff'];
AmpsMedian = [X; medianAmpsDiff'; medianAmpsDevDiff'; medianAmpsNegDiff'; medianAmpsNegDevDiff'];
RTsMedian = [X; medianRTsDiff'; medianRTsDevDiff'; medianRTsNegDiff'; medianRTsNegDevDiff'];
TwoDsMedian = [X; medianTwoDsDiff'; medianTwoDsDevDiff'; medianTwoDsNegDiff'; medianTwoDsNegDevDiff'];
AmpsMin = [X; minAmpsDiff'; minAmpsDevDiff'; minAmpsNegDiff'; minAmpsNegDevDiff'];
RTsMin = [X; minRTsDiff'; minRTsDevDiff'; minRTsNegDiff'; minRTsNegDevDiff'];
TwoDsMin = [X; minTwoDsDiff'; minTwoDsDevDiff'; minTwoDsNegDiff'; minTwoDsNegDevDiff'];
AmpsMax = [X; maxAmpsDiff'; maxAmpsDevDiff'; maxAmpsNegDiff'; maxAmpsNegDevDiff'];
RTsMax = [X; maxRTsDiff'; maxRTsDevDiff'; maxRTsNegDiff'; maxRTsNegDevDiff'];
TwoDsMax = [X; maxTwoDsDiff'; maxTwoDsDevDiff'; maxTwoDsNegDiff'; maxTwoDsNegDevDiff'];
AmpsPrct = [X; prctAmpsDiff'; prctAmpsDevDiff'; prctAmpsNegDiff'; prctAmpsNegDevDiff'];
RTsPrct = [X; prctRTsDiff'; prctRTsDevDiff'; prctRTsNegDiff'; prctRTsNegDevDiff'];
TwoDsPrct = [X; prctTwoDsDiff'; prctTwoDsDevDiff'; prctTwoDsNegDiff'; prctTwoDsNegDevDiff'];

optimData.Amps = optimAmps;
optimData.AmpsNeg = optimAmpsNeg;
optimData.RTs = optimRTs;
optimData.RTsNeg = optimRTsNeg;
optimData.TwoDs = optimTwoDs;
optimData.TwoDsNeg = optimTwoDsNeg;
optimData.SD = SD;

[AmpsMeanBottom, AmpsMedianBottom, AmpsMinBottom, AmpsMaxBottom, AmpsPrctBottom, meanAmpsDiffBottom, meanAmpsDevDiffBottom, meanAmpsNegDiffBottom,...
    meanAmpsNegDevDiffBottom, medianAmpsDiffBottom, medianAmpsDevDiffBottom, medianAmpsNegDiffBottom, medianAmpsNegDevDiffBottom, prctAmpsDiffBottom,...
    prctAmpsDevDiffBottom, prctAmpsNegDiffBottom, prctAmpsNegDevDiffBottom, minAmpsDiffBottom, minAmpsDevDiffBottom, minAmpsNegDiffBottom,...
    minAmpsNegDevDiffBottom, maxAmpsDiffBottom, maxAmpsDevDiffBottom, maxAmpsNegDiffBottom, maxAmpsNegDevDiffBottom, initAmpsBottom,...
    initAmpsNegBottom] = errorAmps(sweepcount, initAmps, initAmpsNeg, prct, prctNeg, Bottom);
[AmpsMeanMid, AmpsMedianMid, AmpsMinMid, AmpsMaxMid, AmpsPrctMid, meanAmpsDiffMid, meanAmpsDevDiffMid, meanAmpsNegDiffMid, meanAmpsNegDevDiffMid,...
    medianAmpsDiffMid, medianAmpsDevDiffMid, medianAmpsNegDiffMid, medianAmpsNegDevDiffMid, prctAmpsDiffMid, prctAmpsDevDiffMid, prctAmpsNegDiffMid,...
    prctAmpsNegDevDiffMid, minAmpsDiffMid, minAmpsDevDiffMid, minAmpsNegDiffMid, minAmpsNegDevDiffMid, maxAmpsDiffMid, maxAmpsDevDiffMid, maxAmpsNegDiffMid,...
    maxAmpsNegDevDiffMid, initAmpsMid, initAmpsNegMid] = errorAmps(sweepcount, initAmps, initAmpsNeg, prct, prctNeg, Mid);
[AmpsMeanTop, AmpsMedianTop, AmpsMinTop, AmpsMaxTop, AmpsPrctTop, meanAmpsDiffTop, meanAmpsDevDiffTop, meanAmpsNegDiffTop, meanAmpsNegDevDiffTop,...
    medianAmpsDiffTop, medianAmpsDevDiffTop, medianAmpsNegDiffTop, medianAmpsNegDevDiffTop, prctAmpsDiffTop, prctAmpsDevDiffTop, prctAmpsNegDiffTop,...
    prctAmpsNegDevDiffTop, minAmpsDiffTop, minAmpsDevDiffTop, minAmpsNegDiffTop, minAmpsNegDevDiffTop, maxAmpsDiffTop, maxAmpsDevDiffTop, maxAmpsNegDiffTop,...
    maxAmpsNegDevDiffTop, initAmpsTop, initAmpsNegTop] = errorAmps(sweepcount, initAmps, initAmpsNeg, prct, prctNeg, Top);





%% Fit curves to data:
warning('off', 'MATLAB:nearlySingularMatrix');

pow = .5;
powWgh = .5;



% Amplitude regular and weighted non-linear fits:
[AmpsErr, AmpsErrWgh, XAmpsErr, XAmpsErrWgh, AmpsErrLSE, AmpsErrWghLSE, XExt] = curvefitMinisSingle(pow, powWgh, X, prctAmpsDiff', nSweeps);
[AmpsDevErr, AmpsDevErrWgh, XAmpsDevErr, XAmpsDevErrWgh, AmpsDevErrLSE, AmpsDevErrWghLSE] = curvefitMinisSingle(pow, powWgh, X, prctAmpsDevDiff', nSweeps);
[AmpsNegErr, AmpsNegErrWgh, XAmpsNegErr, XAmpsNegErrWgh, AmpsNegErrLSE, AmpsNegErrWghLSE] = curvefitMinisSingle(pow, powWgh, X, prctAmpsNegDiff', nSweeps);
[AmpsNegDevErr, AmpsNegDevErrWgh, XAmpsNegDevErr, XAmpsNegDevErrWgh, AmpsNegDevErrLSE, AmpsNegDevErrWghLSE] = curvefitMinisSingle(pow, powWgh, X,...
    prctAmpsNegDevDiff', nSweeps);

% Amplitude regular unconstrained linear fits:
regLineWeights = ones(size(nSweeps));
[AmpsErrRegrUncons, AmpsErrRegrUnconsLSE] = lineMinisUncons(X, prctAmpsDiff', regLineWeights);
[AmpsDevErrRegrUncons, AmpsDevErrRegrUnconsLSE] = lineMinisUncons(X, prctAmpsDevDiff', regLineWeights);
[AmpsNegErrRegrUncons, AmpsNegErrRegrUnconsLSE] = lineMinisUncons(X, prctAmpsNegDiff', regLineWeights);
[AmpsNegDevErrRegrUncons, AmpsNegDevErrRegrUnconsLSE] = lineMinisUncons(X, prctAmpsNegDevDiff', regLineWeights);

% Amplitude weighted unconstrained linear fits:
[AmpsErrWghRegrUncons, AmpsErrWghRegrUnconsLSE] = lineMinisUncons(X, prctAmpsDiff', nSweeps);
[AmpsDevErrWghRegrUncons, AmpsDevErrWghRegrUnconsLSE] = lineMinisUncons(X, prctAmpsDevDiff', nSweeps);
[AmpsNegErrWghRegrUncons, AmpsNegErrWghRegrUnconsLSE] = lineMinisUncons(X, prctAmpsNegDiff', nSweeps);
[AmpsNegDevErrWghRegrUncons, AmpsNegDevErrWghRegrUnconsLSE] = lineMinisUncons(X, prctAmpsNegDevDiff', nSweeps);

% Bottom amplitude regular and weighted non-linear fits:
[AmpsErrBottom, AmpsErrWghBottom, XAmpsErrBottom, XAmpsErrWghBottom, AmpsErrLSEBottom, AmpsErrWghLSEBottom] = curvefitMinisSingle(pow, powWgh, X,...
    prctAmpsDiffBottom', nSweeps);
[AmpsDevErrBottom, AmpsDevErrWghBottom, XAmpsDevErrBottom, XAmpsDevErrWghBottom, AmpsDevErrLSEBottom, AmpsDevErrWghLSEBottom] = curvefitMinisSingle(pow,...
    powWgh, X, prctAmpsDevDiffBottom', nSweeps);
[AmpsNegErrBottom, AmpsNegErrWghBottom, XAmpsNegErrBottom, XAmpsNegErrWghBottom, AmpsNegErrLSEBottom, AmpsNegErrWghLSEBottom] = curvefitMinisSingle(pow,...
    powWgh, X, prctAmpsNegDiffBottom', nSweeps);
[AmpsNegDevErrBottom, AmpsNegDevErrWghBottom, XAmpsNegDevErrBottom, XAmpsNegDevErrWghBottom, AmpsNegDevErrLSEBottom,...
    AmpsNegDevErrWghLSEBottom] = curvefitMinisSingle(pow, powWgh, X, prctAmpsNegDevDiffBottom', nSweeps);

% Bottom amplitude regular unconstrained linear fits:
[AmpsErrRegrUnconsBottom, AmpsErrRegrUnconsLSEBottom] = lineMinisUncons(X, prctAmpsDiffBottom', regLineWeights);
[AmpsDevErrRegrUnconsBottom, AmpsDevErrRegrUnconsLSEBottom] = lineMinisUncons(X, prctAmpsDevDiffBottom', regLineWeights);
[AmpsNegErrRegrUnconsBottom, AmpsNegErrRegrUnconsLSEBottom] = lineMinisUncons(X, prctAmpsNegDiffBottom', regLineWeights);
[AmpsNegDevErrRegrUnconsBottom, AmpsNegDevErrRegrUnconsLSEBottom] = lineMinisUncons(X, prctAmpsNegDevDiffBottom', regLineWeights);

% Bottom amplitude weighted unconstrained linear fits:
[AmpsErrWghRegrUnconsBottom, AmpsErrWghRegrUnconsLSEBottom] = lineMinisUncons(X, prctAmpsDiffBottom', nSweeps);
[AmpsDevErrWghRegrUnconsBottom, AmpsDevErrWghRegrUnconsLSEBottom] = lineMinisUncons(X, prctAmpsDevDiffBottom', nSweeps);
[AmpsNegErrWghRegrUnconsBottom, AmpsNegErrWghRegrUnconsLSEBottom] = lineMinisUncons(X, prctAmpsNegDiffBottom', nSweeps);
[AmpsNegDevErrWghRegrUnconsBottom, AmpsNegDevErrWghRegrUnconsLSEBottom] = lineMinisUncons(X, prctAmpsNegDevDiffBottom', nSweeps);

% Middle amplitude regular and weighted non-linear fits:
[AmpsErrMid, AmpsErrWghMid, XAmpsErrMid, XAmpsErrWghMid, AmpsErrLSEMid, AmpsErrWghLSEMid] = curvefitMinisSingle(pow, powWgh, X, prctAmpsDiffMid', nSweeps);
[AmpsDevErrMid, AmpsDevErrWghMid, XAmpsDevErrMid, XAmpsDevErrWghMid, AmpsDevErrLSEMid, AmpsDevErrWghLSEMid] = curvefitMinisSingle(pow, powWgh, X,...
    prctAmpsDevDiffMid', nSweeps);
[AmpsNegErrMid, AmpsNegErrWghMid, XAmpsNegErrMid, XAmpsNegErrWghMid, AmpsNegErrLSEMid, AmpsNegErrWghLSEMid] = curvefitMinisSingle(pow, powWgh, X,...
    prctAmpsNegDiffMid', nSweeps);
[AmpsNegDevErrMid, AmpsNegDevErrWghMid, XAmpsNegDevErrMid, XAmpsNegDevErrWghMid, AmpsNegDevErrLSEMid, AmpsNegDevErrWghLSEMid] = curvefitMinisSingle(pow,...
    powWgh, X, prctAmpsNegDevDiffMid', nSweeps);

% Middle amplitude regular unconstrained linear fits:
[AmpsErrRegrUnconsMid, AmpsErrRegrUnconsLSEMid] = lineMinisUncons(X, prctAmpsDiffMid', regLineWeights);
[AmpsDevErrRegrUnconsMid, AmpsDevErrRegrUnconsLSEMid] = lineMinisUncons(X, prctAmpsDevDiffMid', regLineWeights);
[AmpsNegErrRegrUnconsMid, AmpsNegErrRegrUnconsLSEMid] = lineMinisUncons(X, prctAmpsNegDiffMid', regLineWeights);
[AmpsNegDevErrRegrUnconsMid, AmpsNegDevErrRegrUnconsLSEMid] = lineMinisUncons(X, prctAmpsNegDevDiffMid', regLineWeights);

% Middle amplitude weighted unconstrained linear fits:
[AmpsErrWghRegrUnconsMid, AmpsErrWghRegrUnconsLSEMid] = lineMinisUncons(X, prctAmpsDiffMid', nSweeps);
[AmpsDevErrWghRegrUnconsMid, AmpsDevErrWghRegrUnconsLSEMid] = lineMinisUncons(X, prctAmpsDevDiffMid', nSweeps);
[AmpsNegErrWghRegrUnconsMid, AmpsNegErrWghRegrUnconsLSEMid] = lineMinisUncons(X, prctAmpsNegDiffMid', nSweeps);
[AmpsNegDevErrWghRegrUnconsMid, AmpsNegDevErrWghRegrUnconsLSEMid] = lineMinisUncons(X, prctAmpsNegDevDiffMid', nSweeps);

% Top amplitude regular and weighted non-linear fits:
[AmpsErrTop, AmpsErrWghTop, XAmpsErrTop, XAmpsErrWghTop, AmpsErrLSETop, AmpsErrWghLSETop] = curvefitMinisSingle(pow, powWgh, X, prctAmpsDiffTop', nSweeps);
[AmpsDevErrTop, AmpsDevErrWghTop, XAmpsDevErrTop, XAmpsDevErrWghTop, AmpsDevErrLSETop, AmpsDevErrWghLSETop] = curvefitMinisSingle(pow, powWgh, X,...
    prctAmpsDevDiffTop', nSweeps);
[AmpsNegErrTop, AmpsNegErrWghTop, XAmpsNegErrTop, XAmpsNegErrWghTop, AmpsNegErrLSETop, AmpsNegErrWghLSETop] = curvefitMinisSingle(pow, powWgh, X,...
    prctAmpsNegDiffTop', nSweeps);
[AmpsNegDevErrTop, AmpsNegDevErrWghTop, XAmpsNegDevErrTop, XAmpsNegDevErrWghTop, AmpsNegDevErrLSETop, AmpsNegDevErrWghLSETop] = curvefitMinisSingle(pow,...
    powWgh, X, prctAmpsNegDevDiffTop', nSweeps);

% Top amplitude regular unconstrained linear fits:
[AmpsErrRegrUnconsTop, AmpsErrRegrUnconsLSETop] = lineMinisUncons(X, prctAmpsDiffTop', regLineWeights);
[AmpsDevErrRegrUnconsTop, AmpsDevErrRegrUnconsLSETop] = lineMinisUncons(X, prctAmpsDevDiffTop', regLineWeights);
[AmpsNegErrRegrUnconsTop, AmpsNegErrRegrUnconsLSETop] = lineMinisUncons(X, prctAmpsNegDiffTop', regLineWeights);
[AmpsNegDevErrRegrUnconsTop, AmpsNegDevErrRegrUnconsLSETop] = lineMinisUncons(X, prctAmpsNegDevDiffTop', regLineWeights);

% Top amplitude weighted unconstrained linear fits:
[AmpsErrWghRegrUnconsTop, AmpsErrWghRegrUnconsLSETop] = lineMinisUncons(X, prctAmpsDiffTop', nSweeps);
[AmpsDevErrWghRegrUnconsTop, AmpsDevErrWghRegrUnconsLSETop] = lineMinisUncons(X, prctAmpsDevDiffTop', nSweeps);
[AmpsNegErrWghRegrUnconsTop, AmpsNegErrWghRegrUnconsLSETop] = lineMinisUncons(X, prctAmpsNegDiffTop', nSweeps);
[AmpsNegDevErrWghRegrUnconsTop, AmpsNegDevErrWghRegrUnconsLSETop] = lineMinisUncons(X, prctAmpsNegDevDiffTop', nSweeps);



% Rise time regular and weighted non-linear fits:
[RTsErr, RTsErrWgh, XRTsErr, XRTsErrWgh, RTsErrLSE, RTsErrWghLSE] = curvefitMinisSingle(pow, powWgh, X, prctRTsDiff', nSweeps);
[RTsDevErr, RTsDevErrWgh, XRTsDevErr, XRTsDevErrWgh, RTsDevErrLSE, RTsDevErrWghLSE] = curvefitMinisSingle(pow, powWgh, X, prctRTsDevDiff', nSweeps);
[RTsNegErr, RTsNegErrWgh, XRTsNegErr, XRTsNegErrWgh, RTsNegErrLSE, RTsNegErrWghLSE] = curvefitMinisSingle(pow, powWgh, X, prctRTsNegDiff', nSweeps);
[RTsNegDevErr, RTsNegDevErrWgh, XRTsNegDevErr, XRTsNegDevErrWgh, RTsNegDevErrLSE, RTsNegDevErrWghLSE] = curvefitMinisSingle(pow, powWgh, X,...
    prctRTsNegDevDiff', nSweeps);

% Rise time regular unconstrained linear fits:
[RTsErrRegrUncons, RTsErrRegrUnconsLSE] = lineMinisUncons(X, prctRTsDiff', regLineWeights);
[RTsDevErrRegrUncons, RTsDevErrRegrUnconsLSE] = lineMinisUncons(X, prctRTsDevDiff', regLineWeights);
[RTsNegErrRegrUncons, RTsNegErrRegrUnconsLSE] = lineMinisUncons(X, prctRTsNegDiff', regLineWeights);
[RTsNegDevErrRegrUncons, RTsNegDevErrRegrUnconsLSE] = lineMinisUncons(X, prctRTsNegDevDiff', regLineWeights);

% Rise time weighted unconstrained linear fits:
[RTsErrWghRegrUncons, RTsErrWghRegrUnconsLSE] = lineMinisUncons(X, prctRTsDiff', nSweeps);
[RTsDevErrWghRegrUncons, RTsDevErrWghRegrUnconsLSE] = lineMinisUncons(X, prctRTsDevDiff', nSweeps);
[RTsNegErrWghRegrUncons, RTsNegErrWghRegrUnconsLSE] = lineMinisUncons(X, prctRTsNegDiff', nSweeps);
[RTsNegDevErrWghRegrUncons, RTsNegDevErrWghRegrUnconsLSE] = lineMinisUncons(X, prctRTsNegDevDiff', nSweeps);



% Combined amplitude and rise time regular and weighted non-linear fits:
[TwoDsErr, TwoDsErrWgh, XTwoDsErr, XTwoDsErrWgh, TwoDsErrLSE, TwoDsErrWghLSE] = curvefitMinisSingle(pow, powWgh, X, prctTwoDsDiff', nSweeps);
[TwoDsDevErr, TwoDsDevErrWgh, XTwoDsDevErr, XTwoDsDevErrWgh, TwoDsDevErrLSE, TwoDsDevErrWghLSE] = curvefitMinisSingle(pow, powWgh, X, prctTwoDsDevDiff', nSweeps);
[TwoDsNegErr, TwoDsNegErrWgh, XTwoDsNegErr, XTwoDsNegErrWgh, TwoDsNegErrLSE, TwoDsNegErrWghLSE] = curvefitMinisSingle(pow, powWgh, X, prctTwoDsNegDiff', nSweeps);
[TwoDsNegDevErr, TwoDsNegDevErrWgh, XTwoDsNegDevErr, XTwoDsNegDevErrWgh, TwoDsNegDevErrLSE, TwoDsNegDevErrWghLSE] = curvefitMinisSingle(pow, powWgh, X,...
    prctTwoDsNegDevDiff', nSweeps);

% Combined amplitude and rise time regular unconstrained linear fits:
[TwoDsErrRegrUncons, TwoDsErrRegrUnconsLSE] = lineMinisUncons(X, prctTwoDsDiff', regLineWeights);
[TwoDsDevErrRegrUncons, TwoDsDevErrRegrUnconsLSE] = lineMinisUncons(X, prctTwoDsDevDiff', regLineWeights);
[TwoDsNegErrRegrUncons, TwoDsNegErrRegrUnconsLSE] = lineMinisUncons(X, prctTwoDsNegDiff', regLineWeights);
[TwoDsNegDevErrRegrUncons, TwoDsNegDevErrRegrUnconsLSE] = lineMinisUncons(X, prctTwoDsNegDevDiff', regLineWeights);

% Combined amplitude and rise time weighted unconstrained linear fits:
[TwoDsErrWghRegrUncons, TwoDsErrWghRegrUnconsLSE] = lineMinisUncons(X, prctTwoDsDiff', nSweeps);
[TwoDsDevErrWghRegrUncons, TwoDsDevErrWghRegrUnconsLSE] = lineMinisUncons(X, prctTwoDsDevDiff', nSweeps);
[TwoDsNegErrWghRegrUncons, TwoDsNegErrWghRegrUnconsLSE] = lineMinisUncons(X, prctTwoDsNegDiff', nSweeps);
[TwoDsNegDevErrWghRegrUncons, TwoDsNegDevErrWghRegrUnconsLSE] = lineMinisUncons(X, prctTwoDsNegDevDiff', nSweeps);

warning('on', 'MATLAB:nearlySingularMatrix');





%% Plot the 50th centile of combined mini EPSP and IPSP data:
if expand
    centFig = plotCentiles(meanAmpsDiff(1), meanRTsDiff(1), medianAmpsDiff(1), medianRTsDiff(1), diffAmps{1}, diffRTs{1}, prctAmpsDiff(1),...
        prctRTsDiff(1), 'EPSP 6-score combined SAD 50th centile');
    centFigNeg = plotCentiles(meanAmpsNegDiff(1), meanRTsNegDiff(1), medianAmpsNegDiff(1), medianRTsNegDiff(1), diffAmpsNeg{1}, diffRTsNeg{1},...
        prctAmpsNegDiff(1), prctRTsNegDiff(1), 'IPSP 6-score combined SAD 50th centile');
else
    centFig = plotCentiles(meanAmpsDiff(1), meanRTsDiff(1), medianAmpsDiff(1), medianRTsDiff(1), diffAmps{1}, diffRTs{1}, prctAmpsDiff(1),...
        prctRTsDiff(1), 'EPSP 3-score combined SAD 50th centile');
    centFigNeg = plotCentiles(meanAmpsNegDiff(1), meanRTsNegDiff(1), medianAmpsNegDiff(1), medianRTsNegDiff(1), diffAmpsNeg{1}, diffRTsNeg{1},...
        prctAmpsNegDiff(1), prctRTsNegDiff(1), 'IPSP 3-score combined SAD 50th centile');
end
cntF = [centFig centFigNeg];





%% Concatenate error bound estimates:
XExtrapol = 1:length(XAmpsErr);
boundAmps = [XExtrapol; XAmpsErr; XAmpsDevErr; XAmpsNegErr; XAmpsNegDevErr; XAmpsErrWgh; XAmpsDevErrWgh; XAmpsNegErrWgh; XAmpsNegDevErrWgh];
AmpsLSE = [AmpsErrLSE, AmpsDevErrLSE, AmpsNegErrLSE, AmpsNegDevErrLSE, AmpsErrWghLSE, AmpsDevErrWghLSE, AmpsNegErrWghLSE, AmpsNegDevErrWghLSE];

boundAmpsBottom = [XExtrapol; XAmpsErrBottom; XAmpsDevErrBottom; XAmpsNegErrBottom; XAmpsNegDevErrBottom; XAmpsErrWghBottom; XAmpsDevErrWghBottom;...
    XAmpsNegErrWghBottom; XAmpsNegDevErrWghBottom];
AmpsLSEBottom = [AmpsErrLSEBottom, AmpsDevErrLSEBottom, AmpsNegErrLSEBottom, AmpsNegDevErrLSEBottom, AmpsErrWghLSEBottom, AmpsDevErrWghLSEBottom,...
    AmpsNegErrWghLSEBottom, AmpsNegDevErrWghLSEBottom];

boundAmpsMid = [XExtrapol; XAmpsErrMid; XAmpsDevErrMid; XAmpsNegErrMid; XAmpsNegDevErrMid; XAmpsErrWghMid; XAmpsDevErrWghMid; XAmpsNegErrWghMid;...
    XAmpsNegDevErrWghMid];
AmpsLSEMid = [AmpsErrLSEMid, AmpsDevErrLSEMid, AmpsNegErrLSEMid, AmpsNegDevErrLSEMid, AmpsErrWghLSEMid, AmpsDevErrWghLSEMid, AmpsNegErrWghLSEMid,...
    AmpsNegDevErrWghLSEMid];

boundAmpsTop = [XExtrapol; XAmpsErrTop; XAmpsDevErrTop; XAmpsNegErrTop; XAmpsNegDevErrTop; XAmpsErrWghTop; XAmpsDevErrWghTop; XAmpsNegErrWghTop;...
    XAmpsNegDevErrWghTop];
AmpsLSETop = [AmpsErrLSETop, AmpsDevErrLSETop, AmpsNegErrLSETop, AmpsNegDevErrLSETop, AmpsErrWghLSETop, AmpsDevErrWghLSETop, AmpsNegErrWghLSETop,...
    AmpsNegDevErrWghLSETop];

boundRTs = [XExtrapol; XRTsErr; XRTsDevErr; XRTsNegErr; XRTsNegDevErr; XRTsErrWgh; XRTsDevErrWgh; XRTsNegErrWgh; XRTsNegDevErrWgh];
RTsLSE = [RTsErrLSE, RTsDevErrLSE, RTsNegErrLSE, RTsNegDevErrLSE, RTsErrWghLSE, RTsDevErrWghLSE, RTsNegErrWghLSE, RTsNegDevErrWghLSE];

boundTwoDs = [XExtrapol; XTwoDsErr; XTwoDsDevErr; XTwoDsNegErr; XTwoDsNegDevErr; XTwoDsErrWgh; XTwoDsDevErrWgh; XTwoDsNegErrWgh; XTwoDsNegDevErrWgh];
TwoDsLSE = [TwoDsErrLSE, TwoDsDevErrLSE, TwoDsNegErrLSE, TwoDsNegDevErrLSE, TwoDsErrWghLSE, TwoDsDevErrWghLSE, TwoDsNegErrWghLSE, TwoDsNegDevErrWghLSE];

boundAmpsLinear = [XExtrapol; AmpsErrRegrUncons; AmpsDevErrRegrUncons; AmpsNegErrRegrUncons; AmpsNegDevErrRegrUncons;...
    AmpsErrWghRegrUncons; AmpsDevErrWghRegrUncons; AmpsNegErrWghRegrUncons; AmpsNegDevErrWghRegrUncons];
AmpsLinearLSE = [AmpsErrRegrUnconsLSE, AmpsDevErrRegrUnconsLSE, AmpsNegErrRegrUnconsLSE, AmpsNegDevErrRegrUnconsLSE,...
    AmpsErrWghRegrUnconsLSE, AmpsDevErrWghRegrUnconsLSE, AmpsNegErrWghRegrUnconsLSE, AmpsNegDevErrWghRegrUnconsLSE];

boundAmpsLinearBottom = [XExtrapol; AmpsErrRegrUnconsBottom; AmpsDevErrRegrUnconsBottom; AmpsNegErrRegrUnconsBottom; AmpsNegDevErrRegrUnconsBottom;...
    AmpsErrWghRegrUnconsBottom; AmpsDevErrWghRegrUnconsBottom; AmpsNegErrWghRegrUnconsBottom; AmpsNegDevErrWghRegrUnconsBottom];
AmpsLinearLSEBottom = [AmpsErrRegrUnconsLSEBottom, AmpsDevErrRegrUnconsLSEBottom, AmpsNegErrRegrUnconsLSEBottom, AmpsNegDevErrRegrUnconsLSEBottom,...
    AmpsErrWghRegrUnconsLSEBottom, AmpsDevErrWghRegrUnconsLSEBottom, AmpsNegErrWghRegrUnconsLSEBottom, AmpsNegDevErrWghRegrUnconsLSEBottom];

boundAmpsLinearMid = [XExtrapol; AmpsErrRegrUnconsMid; AmpsDevErrRegrUnconsMid; AmpsNegErrRegrUnconsMid; AmpsNegDevErrRegrUnconsMid;...
    AmpsErrWghRegrUnconsMid; AmpsDevErrWghRegrUnconsMid; AmpsNegErrWghRegrUnconsMid; AmpsNegDevErrWghRegrUnconsMid];
AmpsLinearLSEMid = [AmpsErrRegrUnconsLSEMid, AmpsDevErrRegrUnconsLSEMid, AmpsNegErrRegrUnconsLSEMid, AmpsNegDevErrRegrUnconsLSEMid,...
    AmpsErrWghRegrUnconsLSEMid, AmpsDevErrWghRegrUnconsLSEMid, AmpsNegErrWghRegrUnconsLSEMid, AmpsNegDevErrWghRegrUnconsLSEMid];

boundAmpsLinearTop = [XExtrapol; AmpsErrRegrUnconsTop; AmpsDevErrRegrUnconsTop; AmpsNegErrRegrUnconsTop; AmpsNegDevErrRegrUnconsTop;...
    AmpsErrWghRegrUnconsTop; AmpsDevErrWghRegrUnconsTop; AmpsNegErrWghRegrUnconsTop; AmpsNegDevErrWghRegrUnconsTop];
AmpsLinearLSETop = [AmpsErrRegrUnconsLSETop, AmpsDevErrRegrUnconsLSETop, AmpsNegErrRegrUnconsLSETop, AmpsNegDevErrRegrUnconsLSETop,...
    AmpsErrWghRegrUnconsLSETop, AmpsDevErrWghRegrUnconsLSETop, AmpsNegErrWghRegrUnconsLSETop, AmpsNegDevErrWghRegrUnconsLSETop];

boundRTsLinear = [XExtrapol; RTsErrRegrUncons; RTsDevErrRegrUncons; RTsNegErrRegrUncons; RTsNegDevErrRegrUncons;...
    RTsErrWghRegrUncons; RTsDevErrWghRegrUncons; RTsNegErrWghRegrUncons; RTsNegDevErrWghRegrUncons];
RTsLinearLSE = [RTsErrRegrUnconsLSE, RTsDevErrRegrUnconsLSE, RTsNegErrRegrUnconsLSE, RTsNegDevErrRegrUnconsLSE,...
    RTsErrWghRegrUnconsLSE, RTsDevErrWghRegrUnconsLSE, RTsNegErrWghRegrUnconsLSE, RTsNegDevErrWghRegrUnconsLSE];

boundTwoDsLinear = [XExtrapol; TwoDsErrRegrUncons; TwoDsDevErrRegrUncons; TwoDsNegErrRegrUncons; TwoDsNegDevErrRegrUncons;...
    TwoDsErrWghRegrUncons; TwoDsDevErrWghRegrUncons; TwoDsNegErrWghRegrUncons; TwoDsNegDevErrWghRegrUncons];
TwoDsLinearLSE = [TwoDsErrRegrUnconsLSE, TwoDsDevErrRegrUnconsLSE, TwoDsNegErrRegrUnconsLSE, TwoDsNegDevErrRegrUnconsLSE,...
    TwoDsErrWghRegrUnconsLSE, TwoDsDevErrWghRegrUnconsLSE, TwoDsNegErrWghRegrUnconsLSE, TwoDsNegDevErrWghRegrUnconsLSE];





%% Plot regular non-linear fits:
button = questdlg('Display the model fitting graphs?','Display Graphs','Yes','No','Yes');
if strcmpi(button, 'Yes')
    
    % Amplitude error bounds:
    f1 = plotErrBounds(X, XExt, AmpsErr, prctAmpsDiff, meanAmpsDiff, medianAmpsDiff, minAmpsDiff, maxAmpsDiff, AmpsErrLSE,...
        'Mini EPSP amplitude SAD bound (regular non-linear fit)', [], expand);
    f2 = plotErrBounds(X, XExt, AmpsDevErr, prctAmpsDevDiff, meanAmpsDevDiff, medianAmpsDevDiff, minAmpsDevDiff, maxAmpsDevDiff, AmpsDevErrLSE,...
        'Mini EPSP amplitude MAD bound (regular non-linear fit)', [], expand);
    f3 = plotErrBounds(X, XExt, AmpsNegErr, prctAmpsNegDiff, meanAmpsNegDiff, medianAmpsNegDiff, minAmpsNegDiff, maxAmpsNegDiff, AmpsNegErrLSE,...
        'Mini IPSP amplitude SAD bound (regular non-linear fit)', [], expand);
    f4 = plotErrBounds(X, XExt, AmpsNegDevErr, prctAmpsNegDevDiff, meanAmpsNegDevDiff, medianAmpsNegDevDiff, minAmpsNegDevDiff, maxAmpsNegDevDiff,...
        AmpsNegDevErrLSE, 'Mini IPSP amplitude MAD bound (regular non-linear fit)', [], expand);
    
    % Bottom amplitude error bounds:
    f1Bottom = plotErrBounds(X, XExt, AmpsErrBottom, prctAmpsDiffBottom, meanAmpsDiffBottom, medianAmpsDiffBottom, minAmpsDiffBottom, maxAmpsDiffBottom,...
        AmpsErrLSEBottom, 'Mini EPSP amplitude top 50% SAD bound (regular non-linear fit)', [], expand);
    f2Bottom = plotErrBounds(X, XExt, AmpsDevErrBottom, prctAmpsDevDiffBottom, meanAmpsDevDiffBottom, medianAmpsDevDiffBottom, minAmpsDevDiffBottom,...
        maxAmpsDevDiffBottom, AmpsDevErrLSEBottom, 'Mini EPSP amplitude top 50% MAD bound (regular non-linear fit)', [], expand);
    f3Bottom = plotErrBounds(X, XExt, AmpsNegErrBottom, prctAmpsNegDiffBottom, meanAmpsNegDiffBottom, medianAmpsNegDiffBottom, minAmpsNegDiffBottom,...
        maxAmpsNegDiffBottom, AmpsNegErrLSEBottom, 'Mini IPSP amplitude top 50% SAD bound (regular non-linear fit)', [], expand);
    f4Bottom = plotErrBounds(X, XExt, AmpsNegDevErrBottom, prctAmpsNegDevDiffBottom, meanAmpsNegDevDiffBottom, medianAmpsNegDevDiffBottom,...
        minAmpsNegDevDiffBottom, maxAmpsNegDevDiffBottom, AmpsNegDevErrLSEBottom, 'Mini IPSP amplitude top 50% MAD bound (regular non-linear fit)', [], expand);
    
    % Middle amplitude error bounds:
    f1Mid = plotErrBounds(X, XExt, AmpsErrMid, prctAmpsDiffMid, meanAmpsDiffMid, medianAmpsDiffMid, minAmpsDiffMid, maxAmpsDiffMid, AmpsErrLSEMid,...
        'Mini EPSP amplitude top 10% SAD bound (regular non-linear fit)', [], expand);
    f2Mid = plotErrBounds(X, XExt, AmpsDevErrMid, prctAmpsDevDiffMid, meanAmpsDevDiffMid, medianAmpsDevDiffMid, minAmpsDevDiffMid, maxAmpsDevDiffMid,...
        AmpsDevErrLSEMid, 'Mini EPSP amplitude top 10% MAD bound (regular non-linear fit)', [], expand);
    f3Mid = plotErrBounds(X, XExt, AmpsNegErrMid, prctAmpsNegDiffMid, meanAmpsNegDiffMid, medianAmpsNegDiffMid, minAmpsNegDiffMid, maxAmpsNegDiffMid,...
        AmpsNegErrLSEMid, 'Mini IPSP amplitude top 10% SAD bound (regular non-linear fit)', [], expand);
    f4Mid = plotErrBounds(X, XExt, AmpsNegDevErrMid, prctAmpsNegDevDiffMid, meanAmpsNegDevDiffMid, medianAmpsNegDevDiffMid, minAmpsNegDevDiffMid,...
        maxAmpsNegDevDiffMid, AmpsNegDevErrLSEMid, 'Mini IPSP amplitude top 10% MAD bound (regular non-linear fit)', [], expand);
    
    % Top amplitude error bounds:
    f1Top = plotErrBounds(X, XExt, AmpsErrTop, prctAmpsDiffTop, meanAmpsDiffTop, medianAmpsDiffTop, minAmpsDiffTop, maxAmpsDiffTop, AmpsErrLSETop,...
        'Mini EPSP amplitude top 2% SAD bound (regular non-linear fit)', [], expand);
    f2Top = plotErrBounds(X, XExt, AmpsDevErrTop, prctAmpsDevDiffTop, meanAmpsDevDiffTop, medianAmpsDevDiffTop, minAmpsDevDiffTop, maxAmpsDevDiffTop,...
        AmpsDevErrLSETop, 'Mini EPSP amplitude top 2% MAD bound (regular non-linear fit)', [], expand);
    f3Top = plotErrBounds(X, XExt, AmpsNegErrTop, prctAmpsNegDiffTop, meanAmpsNegDiffTop, medianAmpsNegDiffTop, minAmpsNegDiffTop, maxAmpsNegDiffTop,...
        AmpsNegErrLSETop, 'Mini IPSP amplitude top 2% SAD bound (regular non-linear fit)', [], expand);
    f4Top = plotErrBounds(X, XExt, AmpsNegDevErrTop, prctAmpsNegDevDiffTop, meanAmpsNegDevDiffTop, medianAmpsNegDevDiffTop, minAmpsNegDevDiffTop,...
        maxAmpsNegDevDiffTop, AmpsNegDevErrLSETop, 'Mini IPSP amplitude top 2% MAD bound (regular non-linear fit)', [], expand);
    
    % Rise time error bounds:
    g1 = plotErrBounds(X, XExt, RTsErr, prctRTsDiff, meanRTsDiff, medianRTsDiff, minRTsDiff, maxRTsDiff, RTsErrLSE,...
        'Mini EPSP rise time SAD bound (regular non-linear fit)', [], expand);
    g2 = plotErrBounds(X, XExt, RTsDevErr, prctRTsDevDiff, meanRTsDevDiff, medianRTsDevDiff, minRTsDevDiff, maxRTsDevDiff, RTsDevErrLSE,...
        'Mini EPSP rise time MAD bound (regular non-linear fit)', [], expand);
    g3 = plotErrBounds(X, XExt, RTsNegErr, prctRTsNegDiff, meanRTsNegDiff, medianRTsNegDiff, minRTsNegDiff, maxRTsNegDiff, RTsNegErrLSE,...
        'Mini IPSP rise time SAD bound (regular non-linear fit)', [], expand);
    g4 = plotErrBounds(X, XExt, RTsNegDevErr, prctRTsNegDevDiff, meanRTsNegDevDiff, medianRTsNegDevDiff, minRTsNegDevDiff, maxRTsNegDevDiff,...
        RTsNegDevErrLSE, 'Mini IPSP rise time MAD bound (regular non-linear fit)', [], expand);
    
    % Combined amplitude and rise time error bounds:
    h1 = plotErrBounds(X, XExt, TwoDsErr, prctTwoDsDiff, meanTwoDsDiff, medianTwoDsDiff, minTwoDsDiff, maxTwoDsDiff, TwoDsErrLSE,...
        'Mini EPSP combined amplitude and rise time SAD bound (regular non-linear fit)', [], expand);
    h2 = plotErrBounds(X, XExt, TwoDsDevErr, prctTwoDsDevDiff, meanTwoDsDevDiff, medianTwoDsDevDiff, minTwoDsDevDiff, maxTwoDsDevDiff, TwoDsDevErrLSE,...
        'Mini EPSP combined amplitude and rise time MAD bound (regular non-linear fit)', [], expand);
    h3 = plotErrBounds(X, XExt, TwoDsNegErr, prctTwoDsNegDiff, meanTwoDsNegDiff, medianTwoDsNegDiff, minTwoDsNegDiff, maxTwoDsNegDiff, TwoDsNegErrLSE,...
        'Mini IPSP combined amplitude and rise time SAD bound (regular non-linear fit)', [], expand);
    h4 = plotErrBounds(X, XExt, TwoDsNegDevErr, prctTwoDsNegDevDiff, meanTwoDsNegDevDiff, medianTwoDsNegDevDiff, minTwoDsNegDevDiff, maxTwoDsNegDevDiff,...
        TwoDsNegDevErrLSE, 'Mini IPSP combined amplitude and rise time MAD bound (regular non-linear fit)', [], expand);
    
    
    
    
    
    %% Plot weighted non-linear fits:
    % Amplitude error bounds:
    f5 = plotErrBoundsWgh(X, XExt, AmpsErrWgh, prctAmpsDiff, meanAmpsDiff, medianAmpsDiff, minAmpsDiff, maxAmpsDiff, AmpsErrWghLSE, nSweeps,...
        'Mini EPSP amplitude SAD bound (confidence-weighted non-linear fit)', [], expand);
    f6 = plotErrBoundsWgh(X, XExt, AmpsDevErrWgh, prctAmpsDevDiff, meanAmpsDevDiff, medianAmpsDevDiff, minAmpsDevDiff, maxAmpsDevDiff, AmpsDevErrWghLSE,...
        nSweeps, 'Mini EPSP amplitude MAD bound (confidence-weighted non-linear fit)', [], expand);
    f7 = plotErrBoundsWgh(X, XExt, AmpsNegErrWgh, prctAmpsNegDiff, meanAmpsNegDiff, medianAmpsNegDiff, minAmpsNegDiff, maxAmpsNegDiff, AmpsNegErrWghLSE,...
        nSweeps, 'Mini IPSP amplitude SAD bound (confidence-weighted non-linear fit)', [], expand);
    f8 = plotErrBoundsWgh(X, XExt, AmpsNegDevErrWgh, prctAmpsNegDevDiff, meanAmpsNegDevDiff, medianAmpsNegDevDiff, minAmpsNegDevDiff, maxAmpsNegDevDiff,...
        AmpsNegDevErrWghLSE, nSweeps, 'Mini IPSP amplitude MAD bound (confidence-weighted non-linear fit)', [], expand);
    F = [f1 f2 f3 f4 f5 f6 f7 f8];
    
    % Bottom amplitude error bounds:
    f5Bottom = plotErrBoundsWgh(X, XExt, AmpsErrWghBottom, prctAmpsDiffBottom, meanAmpsDiffBottom, medianAmpsDiffBottom, minAmpsDiffBottom,...
        maxAmpsDiffBottom, AmpsErrWghLSEBottom, nSweeps, 'Mini EPSP amplitude top 50% SAD bound (confidence-weighted non-linear fit)', [], expand);
    f6Bottom = plotErrBoundsWgh(X, XExt, AmpsDevErrWghBottom, prctAmpsDevDiffBottom, meanAmpsDevDiffBottom, medianAmpsDevDiffBottom, minAmpsDevDiffBottom,...
        maxAmpsDevDiffBottom, AmpsDevErrWghLSEBottom, nSweeps, 'Mini EPSP amplitude top 50% MAD bound (confidence-weighted non-linear fit)', [], expand);
    f7Bottom = plotErrBoundsWgh(X, XExt, AmpsNegErrWghBottom, prctAmpsNegDiffBottom, meanAmpsNegDiffBottom, medianAmpsNegDiffBottom, minAmpsNegDiffBottom,...
        maxAmpsNegDiffBottom, AmpsNegErrWghLSEBottom, nSweeps, 'Mini IPSP amplitude top 50% SAD bound (confidence-weighted non-linear fit)', [], expand);
    f8Bottom = plotErrBoundsWgh(X, XExt, AmpsNegDevErrWghBottom, prctAmpsNegDevDiffBottom, meanAmpsNegDevDiffBottom, medianAmpsNegDevDiffBottom,...
        minAmpsNegDevDiffBottom, maxAmpsNegDevDiffBottom, AmpsNegDevErrWghLSEBottom, nSweeps,...
        'Mini IPSP amplitude top 50% MAD bound (confidence-weighted non-linear fit)', [], expand);
    FBottom = [f1Bottom f2Bottom f3Bottom f4Bottom f5Bottom f6Bottom f7Bottom f8Bottom];
    
    % Middle amplitude error bounds:
    f5Mid = plotErrBoundsWgh(X, XExt, AmpsErrWghMid, prctAmpsDiffMid, meanAmpsDiffMid, medianAmpsDiffMid, minAmpsDiffMid, maxAmpsDiffMid,...
        AmpsErrWghLSEMid, nSweeps, 'Mini EPSP amplitude top 10% SAD bound (confidence-weighted non-linear fit)', [], expand);
    f6Mid = plotErrBoundsWgh(X, XExt, AmpsDevErrWghMid, prctAmpsDevDiffMid, meanAmpsDevDiffMid, medianAmpsDevDiffMid, minAmpsDevDiffMid,...
        maxAmpsDevDiffMid, AmpsDevErrWghLSEMid, nSweeps, 'Mini EPSP amplitude top 10% MAD bound (confidence-weighted non-linear fit)', [], expand);
    f7Mid = plotErrBoundsWgh(X, XExt, AmpsNegErrWghMid, prctAmpsNegDiffMid, meanAmpsNegDiffMid, medianAmpsNegDiffMid, minAmpsNegDiffMid,...
        maxAmpsNegDiffMid, AmpsNegErrWghLSEMid, nSweeps, 'Mini IPSP amplitude top 10% SAD bound (confidence-weighted non-linear fit)', [], expand);
    f8Mid = plotErrBoundsWgh(X, XExt, AmpsNegDevErrWghMid, prctAmpsNegDevDiffMid, meanAmpsNegDevDiffMid, medianAmpsNegDevDiffMid, minAmpsNegDevDiffMid,...
        maxAmpsNegDevDiffMid, AmpsNegDevErrWghLSEMid, nSweeps, 'Mini IPSP amplitude top 10% MAD bound (confidence-weighted non-linear fit)', [], expand);
    FMid = [f1Mid f2Mid f3Mid f4Mid f5Mid f6Mid f7Mid f8Mid];
    
    % Top amplitude error bounds:
    f5Top = plotErrBoundsWgh(X, XExt, AmpsErrWghTop, prctAmpsDiffTop, meanAmpsDiffTop, medianAmpsDiffTop, minAmpsDiffTop, maxAmpsDiffTop,...
        AmpsErrWghLSETop, nSweeps, 'Mini EPSP amplitude top 2% SAD bound (confidence-weighted non-linear fit)', [], expand);
    f6Top = plotErrBoundsWgh(X, XExt, AmpsDevErrWghTop, prctAmpsDevDiffTop, meanAmpsDevDiffTop, medianAmpsDevDiffTop, minAmpsDevDiffTop,...
        maxAmpsDevDiffTop, AmpsDevErrWghLSETop, nSweeps, 'Mini EPSP amplitude top 2% MAD bound (confidence-weighted non-linear fit)', [], expand);
    f7Top = plotErrBoundsWgh(X, XExt, AmpsNegErrWghTop, prctAmpsNegDiffTop, meanAmpsNegDiffTop, medianAmpsNegDiffTop, minAmpsNegDiffTop,...
        maxAmpsNegDiffTop, AmpsNegErrWghLSETop, nSweeps, 'Mini IPSP amplitude top 2% SAD bound (confidence-weighted non-linear fit)', [], expand);
    f8Top = plotErrBoundsWgh(X, XExt, AmpsNegDevErrWghTop, prctAmpsNegDevDiffTop, meanAmpsNegDevDiffTop, medianAmpsNegDevDiffTop, minAmpsNegDevDiffTop,...
        maxAmpsNegDevDiffTop, AmpsNegDevErrWghLSETop, nSweeps, 'Mini IPSP amplitude top 2% MAD bound (confidence-weighted non-linear fit)', [], expand);
    FTop = [f1Top f2Top f3Top f4Top f5Top f6Top f7Top f8Top];
    
    % Rise time error bounds:
    g5 = plotErrBoundsWgh(X, XExt, RTsErrWgh, prctRTsDiff, meanRTsDiff, medianRTsDiff, minRTsDiff, maxRTsDiff, RTsErrWghLSE, nSweeps,...
        'Mini EPSP rise time SAD bound (confidence-weighted non-linear fit)', [], expand);
    g6 = plotErrBoundsWgh(X, XExt, RTsDevErrWgh, prctRTsDevDiff, meanRTsDevDiff, medianRTsDevDiff, minRTsDevDiff, maxRTsDevDiff, RTsDevErrWghLSE,...
        nSweeps, 'Mini EPSP rise time MAD bound (confidence-weighted non-linear fit)', [], expand);
    g7 = plotErrBoundsWgh(X, XExt, RTsNegErrWgh, prctRTsNegDiff, meanRTsNegDiff, medianRTsNegDiff, minRTsNegDiff, maxRTsNegDiff, RTsNegErrWghLSE,...
        nSweeps, 'Mini IPSP rise time SAD bound (confidence-weighted non-linear fit)', [], expand);
    g8 = plotErrBoundsWgh(X, XExt, RTsNegDevErrWgh, prctRTsNegDevDiff, meanRTsNegDevDiff, medianRTsNegDevDiff, minRTsNegDevDiff, maxRTsNegDevDiff,...
        RTsNegDevErrWghLSE, nSweeps, 'Mini IPSP rise time MAD bound (confidence-weighted non-linear fit)', [], expand);
    G = [g1 g2 g3 g4 g5 g6 g7 g8];
    
    % Combined amplitude and rise time error bounds:
    h5 = plotErrBoundsWgh(X, XExt, TwoDsErrWgh, prctTwoDsDiff, meanTwoDsDiff, medianTwoDsDiff, minTwoDsDiff, maxTwoDsDiff, TwoDsErrWghLSE, nSweeps,...
        'Mini EPSP combined amplitude and rise time SAD bound (confidence-weighted non-linear fit)', [], expand);
    h6 = plotErrBoundsWgh(X, XExt, TwoDsDevErrWgh, prctTwoDsDevDiff, meanTwoDsDevDiff, medianTwoDsDevDiff, minTwoDsDevDiff, maxTwoDsDevDiff,...
        TwoDsDevErrWghLSE, nSweeps, 'Mini EPSP combined amplitude and rise time MAD bound (confidence-weighted non-linear fit)', [], expand);
    h7 = plotErrBoundsWgh(X, XExt, TwoDsNegErrWgh, prctTwoDsNegDiff, meanTwoDsNegDiff, medianTwoDsNegDiff, minTwoDsNegDiff, maxTwoDsNegDiff,...
        TwoDsNegErrWghLSE, nSweeps, 'Mini IPSP combined amplitude and rise time SAD bound (confidence-weighted non-linear fit)', [], expand);
    h8 = plotErrBoundsWgh(X, XExt, TwoDsNegDevErrWgh, prctTwoDsNegDevDiff, meanTwoDsNegDevDiff, medianTwoDsNegDevDiff, minTwoDsNegDevDiff,...
        maxTwoDsNegDevDiff, TwoDsNegDevErrWghLSE, nSweeps, 'Mini IPSP combined amplitude and rise time MAD bound (confidence-weighted non-linear fit)', [], expand);
    H = [h1 h2 h3 h4 h5 h6 h7 h8];
    
    
    
    
    
    %% Plot unconstrained regular linear fits:
    % Amplitude error bounds:
    a1 = plotErrBounds(X, [], AmpsErrRegrUncons, prctAmpsDiff, meanAmpsDiff, medianAmpsDiff, minAmpsDiff, maxAmpsDiff, AmpsErrRegrUnconsLSE,...
        'Mini EPSP amplitude SAD bound (regular linear fit)', 'linear', expand);
    a2 = plotErrBounds(X, [], AmpsDevErrRegrUncons, prctAmpsDevDiff, meanAmpsDevDiff, medianAmpsDevDiff, minAmpsDevDiff, maxAmpsDevDiff,...
        AmpsDevErrRegrUnconsLSE, 'Mini EPSP amplitude MAD bound (regular linear fit)', 'linear', expand);
    a3 = plotErrBounds(X, [], AmpsNegErrRegrUncons, prctAmpsNegDiff, meanAmpsNegDiff, medianAmpsNegDiff, minAmpsNegDiff, maxAmpsNegDiff,...
        AmpsNegErrRegrUnconsLSE, 'Mini IPSP amplitude SAD bound (regular linear fit)', 'linear', expand);
    a4 = plotErrBounds(X, [], AmpsNegDevErrRegrUncons, prctAmpsNegDevDiff, meanAmpsNegDevDiff, medianAmpsNegDevDiff, minAmpsNegDevDiff,...
        maxAmpsNegDevDiff, AmpsNegDevErrRegrUnconsLSE, 'Mini IPSP amplitude MAD bound (regular linear fit)', 'linear', expand);
    
    % Bottom amplitude error bounds:
    a1Bottom = plotErrBounds(X, [], AmpsErrRegrUnconsBottom, prctAmpsDiffBottom, meanAmpsDiffBottom, medianAmpsDiffBottom, minAmpsDiffBottom,...
        maxAmpsDiffBottom, AmpsErrRegrUnconsLSEBottom, 'Mini EPSP amplitude top 50% SAD bound (regular linear fit)', 'linear', expand);
    a2Bottom = plotErrBounds(X, [], AmpsDevErrRegrUnconsBottom, prctAmpsDevDiffBottom, meanAmpsDevDiffBottom, medianAmpsDevDiffBottom,...
        minAmpsDevDiffBottom, maxAmpsDevDiffBottom, AmpsDevErrRegrUnconsLSEBottom, 'Mini EPSP amplitude top 50% MAD bound (regular linear fit)', 'linear', expand);
    a3Bottom = plotErrBounds(X, [], AmpsNegErrRegrUnconsBottom, prctAmpsNegDiffBottom, meanAmpsNegDiffBottom, medianAmpsNegDiffBottom,...
        minAmpsNegDiffBottom, maxAmpsNegDiffBottom, AmpsNegErrRegrUnconsLSEBottom, 'Mini IPSP amplitude top 50% SAD bound (regular linear fit)', 'linear', expand);
    a4Bottom = plotErrBounds(X, [], AmpsNegDevErrRegrUnconsBottom, prctAmpsNegDevDiffBottom, meanAmpsNegDevDiffBottom, medianAmpsNegDevDiffBottom,...
        minAmpsNegDevDiffBottom, maxAmpsNegDevDiffBottom, AmpsNegDevErrRegrUnconsLSEBottom, 'Mini IPSP amplitude top 50% MAD bound (regular linear fit)',...
        'linear', expand);
    
    % Middle amplitude error bounds:
    a1Mid = plotErrBounds(X, [], AmpsErrRegrUnconsMid, prctAmpsDiffMid, meanAmpsDiffMid, medianAmpsDiffMid, minAmpsDiffMid, maxAmpsDiffMid,...
        AmpsErrRegrUnconsLSEMid, 'Mini EPSP amplitude top 10% SAD bound (regular linear fit)', 'linear', expand);
    a2Mid = plotErrBounds(X, [], AmpsDevErrRegrUnconsMid, prctAmpsDevDiffMid, meanAmpsDevDiffMid, medianAmpsDevDiffMid, minAmpsDevDiffMid,...
        maxAmpsDevDiffMid, AmpsDevErrRegrUnconsLSEMid, 'Mini EPSP amplitude top 10% MAD bound (regular linear fit)', 'linear', expand);
    a3Mid = plotErrBounds(X, [], AmpsNegErrRegrUnconsMid, prctAmpsNegDiffMid, meanAmpsNegDiffMid, medianAmpsNegDiffMid, minAmpsNegDiffMid,...
        maxAmpsNegDiffMid, AmpsNegErrRegrUnconsLSEMid, 'Mini IPSP amplitude top 10% SAD bound (regular linear fit)', 'linear', expand);
    a4Mid = plotErrBounds(X, [], AmpsNegDevErrRegrUnconsMid, prctAmpsNegDevDiffMid, meanAmpsNegDevDiffMid, medianAmpsNegDevDiffMid, minAmpsNegDevDiffMid,...
        maxAmpsNegDevDiffMid, AmpsNegDevErrRegrUnconsLSEMid, 'Mini IPSP amplitude top 10% MAD bound (regular linear fit)', 'linear', expand);
    
    % Top amplitude error bounds:
    a1Top = plotErrBounds(X, [], AmpsErrRegrUnconsTop, prctAmpsDiffTop, meanAmpsDiffTop, medianAmpsDiffTop, minAmpsDiffTop, maxAmpsDiffTop,...
        AmpsErrRegrUnconsLSETop, 'Mini EPSP amplitude top 2% SAD bound (regular linear fit)', 'linear', expand);
    a2Top = plotErrBounds(X, [], AmpsDevErrRegrUnconsTop, prctAmpsDevDiffTop, meanAmpsDevDiffTop, medianAmpsDevDiffTop, minAmpsDevDiffTop,...
        maxAmpsDevDiffTop, AmpsDevErrRegrUnconsLSETop, 'Mini EPSP amplitude top 2% MAD bound (regular linear fit)', 'linear', expand);
    a3Top = plotErrBounds(X, [], AmpsNegErrRegrUnconsTop, prctAmpsNegDiffTop, meanAmpsNegDiffTop, medianAmpsNegDiffTop, minAmpsNegDiffTop,...
        maxAmpsNegDiffTop, AmpsNegErrRegrUnconsLSETop, 'Mini IPSP amplitude top 2% SAD bound (regular linear fit)', 'linear', expand);
    a4Top = plotErrBounds(X, [], AmpsNegDevErrRegrUnconsTop, prctAmpsNegDevDiffTop, meanAmpsNegDevDiffTop, medianAmpsNegDevDiffTop, minAmpsNegDevDiffTop,...
        maxAmpsNegDevDiffTop, AmpsNegDevErrRegrUnconsLSETop, 'Mini IPSP amplitude top 2% MAD bound (regular linear fit)', 'linear', expand);
    
    % Rise time error bounds:
    b1 = plotErrBounds(X, [], RTsErrRegrUncons, prctRTsDiff, meanRTsDiff, medianRTsDiff, minRTsDiff, maxRTsDiff, RTsErrRegrUnconsLSE,...
        'Mini EPSP rise time SAD bound (regular linear fit)', 'linear', expand);
    b2 = plotErrBounds(X, [], RTsDevErrRegrUncons, prctRTsDevDiff, meanRTsDevDiff, medianRTsDevDiff, minRTsDevDiff, maxRTsDevDiff,...
        RTsDevErrRegrUnconsLSE, 'Mini EPSP rise time MAD bound (regular linear fit)', 'linear', expand);
    b3 = plotErrBounds(X, [], RTsNegErrRegrUncons, prctRTsNegDiff, meanRTsNegDiff, medianRTsNegDiff, minRTsNegDiff, maxRTsNegDiff,...
        RTsNegErrRegrUnconsLSE, 'Mini IPSP rise time SAD bound (regular linear fit)', 'linear', expand);
    b4 = plotErrBounds(X, [], RTsNegDevErrRegrUncons, prctRTsNegDevDiff, meanRTsNegDevDiff, medianRTsNegDevDiff, minRTsNegDevDiff, maxRTsNegDevDiff,...
        RTsNegDevErrRegrUnconsLSE, 'Mini IPSP rise time MAD bound (regular linear fit)', 'linear', expand);
    
    % Combined amplitude and rise time error bounds:
    c1 = plotErrBounds(X, [], TwoDsErrRegrUncons, prctTwoDsDiff, meanTwoDsDiff, medianTwoDsDiff, minTwoDsDiff, maxTwoDsDiff, TwoDsErrRegrUnconsLSE,...
        'Mini EPSP combined amplitude and rise time SAD bound (regular linear fit)', 'linear', expand);
    c2 = plotErrBounds(X, [], TwoDsDevErrRegrUncons, prctTwoDsDevDiff, meanTwoDsDevDiff, medianTwoDsDevDiff, minTwoDsDevDiff, maxTwoDsDevDiff,...
        TwoDsDevErrRegrUnconsLSE, 'Mini EPSP combined amplitude and rise time MAD bound (regular linear fit)', 'linear', expand);
    c3 = plotErrBounds(X, [], TwoDsNegErrRegrUncons, prctTwoDsNegDiff, meanTwoDsNegDiff, medianTwoDsNegDiff, minTwoDsNegDiff, maxTwoDsNegDiff,...
        TwoDsNegErrRegrUnconsLSE, 'Mini IPSP combined amplitude and rise time SAD bound (regular linear fit)', 'linear', expand);
    c4 = plotErrBounds(X, [], TwoDsNegDevErrRegrUncons, prctTwoDsNegDevDiff, meanTwoDsNegDevDiff, medianTwoDsNegDevDiff, minTwoDsNegDevDiff,...
        maxTwoDsNegDevDiff, TwoDsNegDevErrRegrUnconsLSE, 'Mini IPSP combined amplitude and rise time MAD bound (regular linear fit)', 'linear', expand);
    
    
    
    
    
    %% Plot unconstrained weighted linear fits:
    % Amplitude error bounds:
    a5 = plotErrBoundsWgh(X, [], AmpsErrWghRegrUncons, prctAmpsDiff, meanAmpsDiff, medianAmpsDiff, minAmpsDiff, maxAmpsDiff, AmpsErrWghRegrUnconsLSE,...
        nSweeps, 'Mini EPSP amplitude SAD bound (weighted linear fit)', 'linear', expand);
    a6 = plotErrBoundsWgh(X, [], AmpsDevErrWghRegrUncons, prctAmpsDevDiff, meanAmpsDevDiff, medianAmpsDevDiff, minAmpsDevDiff, maxAmpsDevDiff,...
        AmpsDevErrWghRegrUnconsLSE, nSweeps, 'Mini EPSP amplitude MAD bound (weighted linear fit)', 'linear', expand);
    a7 = plotErrBoundsWgh(X, [], AmpsNegErrWghRegrUncons, prctAmpsNegDiff, meanAmpsNegDiff, medianAmpsNegDiff, minAmpsNegDiff, maxAmpsNegDiff,...
        AmpsNegErrWghRegrUnconsLSE, nSweeps, 'Mini IPSP amplitude SAD bound (weighted linear fit)', 'linear', expand);
    a8 = plotErrBoundsWgh(X, [], AmpsNegDevErrWghRegrUncons, prctAmpsNegDevDiff, meanAmpsNegDevDiff, medianAmpsNegDevDiff, minAmpsNegDevDiff,...
        maxAmpsNegDevDiff, AmpsNegDevErrWghRegrUnconsLSE, nSweeps, 'Mini IPSP amplitude MAD bound (weighted linear fit)', 'linear', expand);
    A = [a1 a2 a3 a4 a5 a6 a7 a8];
    
    % Bottom amplitude error bounds:
    a5Bottom = plotErrBoundsWgh(X, [], AmpsErrWghRegrUnconsBottom, prctAmpsDiffBottom, meanAmpsDiffBottom, medianAmpsDiffBottom, minAmpsDiffBottom,...
        maxAmpsDiffBottom, AmpsErrWghRegrUnconsLSEBottom, nSweeps, 'Mini EPSP amplitude top 50% SAD bound (weighted linear fit)', 'linear', expand);
    a6Bottom = plotErrBoundsWgh(X, [], AmpsDevErrWghRegrUnconsBottom, prctAmpsDevDiffBottom, meanAmpsDevDiffBottom, medianAmpsDevDiffBottom,...
        minAmpsDevDiffBottom, maxAmpsDevDiffBottom, AmpsDevErrWghRegrUnconsLSEBottom, nSweeps,...
        'Mini EPSP amplitude top 50% MAD bound (weighted linear fit)', 'linear', expand);
    a7Bottom = plotErrBoundsWgh(X, [], AmpsNegErrWghRegrUnconsBottom, prctAmpsNegDiffBottom, meanAmpsNegDiffBottom, medianAmpsNegDiffBottom,...
        minAmpsNegDiffBottom, maxAmpsNegDiffBottom, AmpsNegErrWghRegrUnconsLSEBottom, nSweeps,...
        'Mini IPSP amplitude top 50% SAD bound (weighted linear fit)', 'linear', expand);
    a8Bottom = plotErrBoundsWgh(X, [], AmpsNegDevErrWghRegrUnconsBottom, prctAmpsNegDevDiffBottom, meanAmpsNegDevDiffBottom, medianAmpsNegDevDiffBottom,...
        minAmpsNegDevDiffBottom, maxAmpsNegDevDiffBottom, AmpsNegDevErrWghRegrUnconsLSEBottom, nSweeps,...
        'Mini IPSP amplitude top 50% MAD bound (weighted linear fit)', 'linear', expand);
    ABottom = [a1Bottom a2Bottom a3Bottom a4Bottom a5Bottom a6Bottom a7Bottom a8Bottom];
    
    % Middle amplitude error bounds:
    a5Mid = plotErrBoundsWgh(X, [], AmpsErrWghRegrUnconsMid, prctAmpsDiffMid, meanAmpsDiffMid, medianAmpsDiffMid, minAmpsDiffMid, maxAmpsDiffMid,...
        AmpsErrWghRegrUnconsLSEMid, nSweeps, 'Mini EPSP amplitude top 10% SAD bound (weighted linear fit)', 'linear', expand);
    a6Mid = plotErrBoundsWgh(X, [], AmpsDevErrWghRegrUnconsMid, prctAmpsDevDiffMid, meanAmpsDevDiffMid, medianAmpsDevDiffMid, minAmpsDevDiffMid,...
        maxAmpsDevDiffMid, AmpsDevErrWghRegrUnconsLSEMid, nSweeps, 'Mini EPSP amplitude top 10% MAD bound (weighted linear fit)', 'linear', expand);
    a7Mid = plotErrBoundsWgh(X, [], AmpsNegErrWghRegrUnconsMid, prctAmpsNegDiffMid, meanAmpsNegDiffMid, medianAmpsNegDiffMid, minAmpsNegDiffMid,...
        maxAmpsNegDiffMid, AmpsNegErrWghRegrUnconsLSEMid, nSweeps, 'Mini IPSP amplitude top 10% SAD bound (weighted linear fit)', 'linear', expand);
    a8Mid = plotErrBoundsWgh(X, [], AmpsNegDevErrWghRegrUnconsMid, prctAmpsNegDevDiffMid, meanAmpsNegDevDiffMid, medianAmpsNegDevDiffMid,...
        minAmpsNegDevDiffMid, maxAmpsNegDevDiffMid, AmpsNegDevErrWghRegrUnconsLSEMid, nSweeps,...
        'Mini IPSP amplitude top 10% MAD bound (weighted linear fit)', 'linear', expand);
    AMid = [a1Mid a2Mid a3Mid a4Mid a5Mid a6Mid a7Mid a8Mid];
    
    % Top amplitude error bounds:
    a5Top = plotErrBoundsWgh(X, [], AmpsErrWghRegrUnconsTop, prctAmpsDiffTop, meanAmpsDiffTop, medianAmpsDiffTop, minAmpsDiffTop, maxAmpsDiffTop,...
        AmpsErrWghRegrUnconsLSETop, nSweeps, 'Mini EPSP amplitude top 2% SAD bound (weighted linear fit)', 'linear', expand);
    a6Top = plotErrBoundsWgh(X, [], AmpsDevErrWghRegrUnconsTop, prctAmpsDevDiffTop, meanAmpsDevDiffTop, medianAmpsDevDiffTop, minAmpsDevDiffTop,...
        maxAmpsDevDiffTop, AmpsDevErrWghRegrUnconsLSETop, nSweeps, 'Mini EPSP amplitude top 2% MAD bound (weighted linear fit)', 'linear', expand);
    a7Top = plotErrBoundsWgh(X, [], AmpsNegErrWghRegrUnconsTop, prctAmpsNegDiffTop, meanAmpsNegDiffTop, medianAmpsNegDiffTop, minAmpsNegDiffTop,...
        maxAmpsNegDiffTop, AmpsNegErrWghRegrUnconsLSETop, nSweeps, 'Mini IPSP amplitude top 2% SAD bound (weighted linear fit)', 'linear', expand);
    a8Top = plotErrBoundsWgh(X, [], AmpsNegDevErrWghRegrUnconsTop, prctAmpsNegDevDiffTop, meanAmpsNegDevDiffTop, medianAmpsNegDevDiffTop,...
        minAmpsNegDevDiffTop, maxAmpsNegDevDiffTop, AmpsNegDevErrWghRegrUnconsLSETop, nSweeps,...
        'Mini IPSP amplitude top 2% MAD bound (weighted linear fit)', 'linear', expand);
    ATop = [a1Top a2Top a3Top a4Top a5Top a6Top a7Top a8Top];
    
    % Rise time error bounds:
    b5 = plotErrBoundsWgh(X, [], RTsErrWghRegrUncons, prctRTsDiff, meanRTsDiff, medianRTsDiff, minRTsDiff, maxRTsDiff, RTsErrWghRegrUnconsLSE, nSweeps,...
        'Mini EPSP rise time SAD bound (weighted linear fit)', 'linear', expand);
    b6 = plotErrBoundsWgh(X, [], RTsDevErrWghRegrUncons, prctRTsDevDiff, meanRTsDevDiff, medianRTsDevDiff, minRTsDevDiff, maxRTsDevDiff,...
        RTsDevErrWghRegrUnconsLSE, nSweeps, 'Mini EPSP rise time MAD bound (weighted linear fit)', 'linear', expand);
    b7 = plotErrBoundsWgh(X, [], RTsNegErrWghRegrUncons, prctRTsNegDiff, meanRTsNegDiff, medianRTsNegDiff, minRTsNegDiff, maxRTsNegDiff,...
        RTsNegErrWghRegrUnconsLSE, nSweeps, 'Mini IPSP rise time SAD bound (weighted linear fit)', 'linear', expand);
    b8 = plotErrBoundsWgh(X, [], RTsNegDevErrWghRegrUncons, prctRTsNegDevDiff, meanRTsNegDevDiff, medianRTsNegDevDiff, minRTsNegDevDiff,...
        maxRTsNegDevDiff, RTsNegDevErrWghRegrUnconsLSE, nSweeps, 'Mini IPSP rise time MAD bound (weighted linear fit)', 'linear', expand);
    B = [b1 b2 b3 b4 b5 b6 b7 b8];
    
    % Combined amplitude and rise time error bounds:
    c5 = plotErrBoundsWgh(X, [], TwoDsErrWghRegrUncons, prctTwoDsDiff, meanTwoDsDiff, medianTwoDsDiff, minTwoDsDiff, maxTwoDsDiff,...
        TwoDsErrWghRegrUnconsLSE, nSweeps, 'Mini EPSP combined amplitude and rise time SAD bound (weighted linear fit)', 'linear', expand);
    c6 = plotErrBoundsWgh(X, [], TwoDsDevErrWghRegrUncons, prctTwoDsDevDiff, meanTwoDsDevDiff, medianTwoDsDevDiff, minTwoDsDevDiff, maxTwoDsDevDiff,...
        TwoDsDevErrWghRegrUnconsLSE, nSweeps, 'Mini EPSP combined amplitude and rise time MAD bound (weighted linear fit)', 'linear', expand);
    c7 = plotErrBoundsWgh(X, [], TwoDsNegErrWghRegrUncons, prctTwoDsNegDiff, meanTwoDsNegDiff, medianTwoDsNegDiff, minTwoDsNegDiff, maxTwoDsNegDiff,...
        TwoDsNegErrWghRegrUnconsLSE, nSweeps, 'Mini IPSP combined amplitude and rise time SAD bound (weighted linear fit)', 'linear', expand);
    c8 = plotErrBoundsWgh(X, [], TwoDsNegDevErrWghRegrUncons, prctTwoDsNegDevDiff, meanTwoDsNegDevDiff, medianTwoDsNegDevDiff, minTwoDsNegDevDiff,...
        maxTwoDsNegDevDiff, TwoDsNegDevErrWghRegrUnconsLSE, nSweeps, 'Mini IPSP combined amplitude and rise time MAD bound (weighted linear fit)', 'linear', expand);
    C = [c1 c2 c3 c4 c5 c6 c7 c8];
else
    F = [];
    FBottom = [];
    FMid = [];
    FTop = [];
    G = [];
    H = [];
    A = [];
    ABottom = [];
    AMid = [];
    ATop = [];
    B = [];
    C = [];
end





%% Estimate differences between scaled different size files:
estCount = 0;
UFest = [];
button = questdlg('Estimate error bounds for particular size files?','Estimate Error Bounds','Yes','No','Yes');
if strcmpi(button, 'Yes')
    cont = true;
    while cont
        estCountPrev = estCount;
        proceed = false;
        while ~proceed
            filtfs = inputdlg({'Enter the target file size in sweeps:','Enter the noise file size in sweeps:'},'Choose file lengths:',2,{'10','5'});
            if ~isempty(filtfs)
                proceed = true;
                target = str2double(filtfs{1});
                noise = str2double(filtfs{2});
                if target < 1 || noise < 1
                    uiwait(msgbox('The file length cannot be shorter than a single sweep. Please provide only realistic file lengths.', 'modal'));
                    proceed = false;
                end
            else
                proceed = true;
                cont = false;
            end
        end
        if ~isempty(filtfs)
            [estimated, diffAmps, diffAmpsDev, diffAmpsNeg, diffAmpsNegDev, diffRTs, diffRTsDev, diffRTsNeg, diffRTsNegDev, diffTwoDs, diffTwoDsDev,...
                diffTwoDsNeg, diffTwoDsNegDev, minAmps, minAmpsDev, minAmpsNeg, minAmpsNegDev, minRTs, minRTsDev, minRTsNeg, minRTsNegDev, minTwoDs,...
                minTwoDsDev, minTwoDsNeg, minTwoDsNegDev, maxAmps, maxAmpsDev, maxAmpsNeg, maxAmpsNegDev, maxRTs, maxRTsDev, maxRTsNeg, maxRTsNegDev,...
                maxTwoDs, maxTwoDsDev, maxTwoDsNeg, maxTwoDsNegDev, lengthRatio] = compareUneqSweeps(sweepcount, filtfs, initAmps, initAmpsNeg,...
                initRTs, initRTsNeg, inittwoDs, inittwoDsNeg, prct, prctNeg);
            if estimated
                measuredUF = [diffAmps; diffAmpsDev; diffAmpsNeg; diffAmpsNegDev; diffRTs; diffRTsDev; diffRTsNeg; diffRTsNegDev; diffTwoDs;...
                    diffTwoDsDev; diffTwoDsNeg; diffTwoDsNegDev; minAmps; minAmpsDev; minAmpsNeg; minAmpsNegDev; minRTs; minRTsDev; minRTsNeg;...
                    minRTsNegDev; minTwoDs; minTwoDsDev; minTwoDsNeg; minTwoDsNegDev; maxAmps; maxAmpsDev; maxAmpsNeg; maxAmpsNegDev; maxRTs;...
                    maxRTsDev; maxRTsNeg; maxRTsNegDev; maxTwoDs; maxTwoDsDev; maxTwoDsNeg; maxTwoDsNegDev];
                [diffAmpsBottom, diffAmpsDevBottom, diffAmpsNegBottom, diffAmpsNegDevBottom, minAmpsBottom, minAmpsDevBottom, minAmpsNegBottom,...
                    minAmpsNegDevBottom, maxAmpsBottom, maxAmpsDevBottom, maxAmpsNegBottom, maxAmpsNegDevBottom] = compareSweepsAmps(sweepcount, filtfs,...
                    initAmpsBottom, initAmpsNegBottom, prct, prctNeg);
                measuredUFBottom = [diffAmpsBottom; diffAmpsDevBottom; diffAmpsNegBottom; diffAmpsNegDevBottom; minAmpsBottom; minAmpsDevBottom;...
                    minAmpsNegBottom; minAmpsNegDevBottom; maxAmpsBottom; maxAmpsDevBottom; maxAmpsNegBottom; maxAmpsNegDevBottom];
                [diffAmpsMid, diffAmpsDevMid, diffAmpsNegMid, diffAmpsNegDevMid, minAmpsMid, minAmpsDevMid, minAmpsNegMid, minAmpsNegDevMid,...
                    maxAmpsMid, maxAmpsDevMid, maxAmpsNegMid, maxAmpsNegDevMid] = compareSweepsAmps(sweepcount, filtfs, initAmpsMid, initAmpsNegMid,...
                    prct, prctNeg);
                measuredUFMid = [diffAmpsMid; diffAmpsDevMid; diffAmpsNegMid; diffAmpsNegDevMid; minAmpsMid; minAmpsDevMid; minAmpsNegMid;...
                    minAmpsNegDevMid; maxAmpsMid; maxAmpsDevMid; maxAmpsNegMid; maxAmpsNegDevMid];
                [diffAmpsTop, diffAmpsDevTop, diffAmpsNegTop, diffAmpsNegDevTop, minAmpsTop, minAmpsDevTop, minAmpsNegTop, minAmpsNegDevTop,...
                    maxAmpsTop, maxAmpsDevTop, maxAmpsNegTop, maxAmpsNegDevTop] = compareSweepsAmps(sweepcount, filtfs, initAmpsTop, initAmpsNegTop,...
                    prct, prctNeg);
                measuredUFTop = [diffAmpsTop; diffAmpsDevTop; diffAmpsNegTop; diffAmpsNegDevTop; minAmpsTop; minAmpsDevTop; minAmpsNegTop;...
                    minAmpsNegDevTop; maxAmpsTop; maxAmpsDevTop; maxAmpsNegTop; maxAmpsNegDevTop];
            else
                measuredUF = NaN(36,1);
                measuredUFBottom = NaN(12,1);
                measuredUFMid = NaN(12,1);
                measuredUFTop = NaN(12,1);
            end
            predictedUF = predictComparison(target, noise, pow, lengthRatio, XAmpsErr(1), XAmpsDevErr(1), XAmpsNegErr(1), XAmpsNegDevErr(1),...
                XRTsErr(1), XRTsDevErr(1), XRTsNegErr(1), XRTsNegDevErr(1), XTwoDsErr(1), XTwoDsDevErr(1), XTwoDsNegErr(1), XTwoDsNegDevErr(1));
            predictedUFBottom = predictComparisonAmps(target, noise, pow, lengthRatio, XAmpsErrBottom(1), XAmpsDevErrBottom(1), XAmpsNegErrBottom(1),...
                XAmpsNegDevErrBottom(1));
            predictedUFMid = predictComparisonAmps(target, noise, pow, lengthRatio, XAmpsErrMid(1), XAmpsDevErrMid(1), XAmpsNegErrMid(1),...
                XAmpsNegDevErrMid(1));
            predictedUFTop = predictComparisonAmps(target, noise, pow, lengthRatio, XAmpsErrTop(1), XAmpsDevErrTop(1), XAmpsNegErrTop(1),...
                XAmpsNegDevErrTop(1));
            [predictedUFWgh] = predictComparison(target, noise, powWgh, lengthRatio, XAmpsErrWgh(1), XAmpsDevErrWgh(1), XAmpsNegErrWgh(1),...
                XAmpsNegDevErrWgh(1), XRTsErrWgh(1), XRTsDevErrWgh(1), XRTsNegErrWgh(1), XRTsNegDevErrWgh(1), XTwoDsErrWgh(1), XTwoDsDevErrWgh(1),...
                XTwoDsNegErrWgh(1), XTwoDsNegDevErrWgh(1));
            predictedUFWghBottom = predictComparisonAmps(target, noise, pow, lengthRatio, XAmpsErrWghBottom(1), XAmpsDevErrWghBottom(1), XAmpsNegErrWghBottom(1),...
                XAmpsNegDevErrWghBottom(1));
            predictedUFWghMid = predictComparisonAmps(target, noise, pow, lengthRatio, XAmpsErrWghMid(1), XAmpsDevErrWghMid(1), XAmpsNegErrWghMid(1),...
                XAmpsNegDevErrWghMid(1));
            predictedUFWghTop = predictComparisonAmps(target, noise, pow, lengthRatio, XAmpsErrWghTop(1), XAmpsDevErrWghTop(1), XAmpsNegErrWghTop(1),...
                XAmpsNegDevErrWghTop(1));
            UF = [target noise];
            estCount = estCount + 1;
        else
            estimated = false;
        end
        if estCount > estCountPrev
            field = struct('estimated', estimated, 'measuredUF', measuredUF, 'predictedUF', predictedUF, 'predictedUFWgh', predictedUFWgh,...
                'measuredUFBottom', measuredUFBottom, 'predictedUFBottom', predictedUFBottom, 'predictedUFWghBottom', predictedUFWghBottom,...
                'measuredUFMid', measuredUFMid, 'predictedUFMid', predictedUFMid, 'predictedUFWghMid', predictedUFWghMid,...
                'measuredUFTop', measuredUFTop, 'predictedUFTop', predictedUFTop, 'predictedUFWghTop', predictedUFWghTop, 'UF', UF);
            UFest(estCount).estimate = field; %#ok<AGROW>
        end
        if proceed && cont
            button = questdlg('Do you wish to continue estimating?','Estimate Error Bounds','Yes','No','No');
            if strcmpi(button, 'No')
                cont = false;
            end
        end
    end
end
end





%% Local functions:
function [meanCurve, meanCurveWeighted, meanFit, meanFitWeighted, LSE, LSEwgh, xev] = curvefitMinisSingle(pow, powWgh, X, data, weights)

xev = linspace(X(1),X(end),500)';

%A:
a1 = repmat((1:round(data(1))+100)', 1, length(X));
if ~isempty(pow)
    a = (1:length(X)).^pow;
    a = repmat(a, round(data(1))+100, 1);
    a = a.*a1;
    dataRepA = repmat(data,size(a,1),1);
    % Find a curve to fit the entire data vector:
    squareRes = abs(a - dataRepA).^2;
    leastSquareErr = sum(squareRes,2);
    [LSE, iLeastSquareErr] = min(leastSquareErr);
    meanFit = a(iLeastSquareErr,:);
    % Interpolate:
    try
        P = vander(X')\meanFit';
        meanCurve = polyval(P,xev)';
        error(lastwarn);
    catch %#ok<CTCH>
        meanCurve = spline(X, meanFit, xev)';
    end
    % Extrapolate:
    meanFitExt = meanFit(1)*(length(meanFit) + 1: 6*length(meanFit)).^.5;
    meanFit = [meanFit meanFitExt];
else
    meanCurve = [];
    meanFit = [];
    LSE = [];
end

%B:
if ~isempty(powWgh)
    b = (1:length(X)).^powWgh;
    b = repmat(b, round(data(1))+100, 1);
    b = b.*a1;
    dataRepB = repmat(data,size(b,1),1);
    % Find a curve to fit the entire Weighted data vector:
    squareRes = abs(b - dataRepB).^2;
    weights = repmat(weights, size(b,1), 1);
    squareResW = squareRes.*weights;
    leastSquareErr = sum(squareResW,2);
    [~, iLeastSquareErrWeighted] = min(leastSquareErr);
    LSEwgh = sum(squareRes(iLeastSquareErrWeighted,:));
    meanFitWeighted = b(iLeastSquareErrWeighted,:);
    % Interpolate:
    try
        P = vander(X')\meanFitWeighted';
        meanCurveWeighted = polyval(P,xev)';
        error(lastwarn);
    catch %#ok<CTCH>
        meanCurveWeighted = spline(X, meanFitWeighted, xev)';
    end
    % Extrapolate:
    meanFitWeightedExt = meanFitWeighted(1)*(length(meanFitWeighted) + 1: 6*length(meanFitWeighted)).^.5;
    meanFitWeighted = [meanFitWeighted meanFitWeightedExt];
else
    meanCurveWeighted = [];
    meanFitWeighted = [];
    LSEwgh = [];
end

xev = xev';
end



function f = plotErrBounds(X, XExt, Err, centiles, means, medians, mins, maxs, LSE, nameStr, type, expand)

f = figure('position', [50 50 1200 600]);
if strcmpi(type,'linear')
    p1 = plot(X, Err(1:length(X)), 'r');
    hold on
    p2 = plot(X, 0.9*Err(1:length(X)), 'r--');
    plot(X, 1.1*Err(1:length(X)), 'r--');
    plot(X, 0.8*Err(1:length(X)), 'r--');
    plot(X, 1.2*Err(1:length(X)), 'r--');
    plot(X, 0.7*Err(1:length(X)), 'r--');
    plot(X, 1.3*Err(1:length(X)), 'r--');
else
    p1 = plot(XExt, Err, 'r');
    hold on
    p2 = plot(XExt, [0.9*Err; 1.1*Err], 'r--');
    plot(XExt, [0.8*Err; 1.2*Err], 'r--');
    plot(XExt, [0.7*Err; 1.3*Err], 'r--');
end
p3 = plot(X, [mins'; maxs'], 'o', 'MarkerEdgeColor', 'y', 'MarkerFaceColor', 'y');
p4 = plot(X, medians', 's');
p5 = plot(X, means', '^');
p6 = plot(X, centiles', 'o');

set(f, 'Name', nameStr);
title(nameStr);
xlabel('Number of sweeps');
ylabel('Event count');
if expand
    legend([p1 p2(1) p6 p5 p4 p3(1)], 'Model fit to centiles', '10%,20%,30% deviation', '6-score combined 50th centile', 'Mean', 'Median', 'Range',...
        'Location', 'NorthWest');
else
    legend([p1 p2(1) p6 p5 p4 p3(1)], 'Model fit to centiles', '10%,20%,30% deviation', '3-score combined 50th centile', 'Mean', 'Median', 'Range',...
        'Location', 'NorthWest');
end

XLim = xlim;
xpos = XLim(1) + (XLim(2)-XLim(1))*.01;
YLim = ylim;
ypos = YLim(2) - (YLim(2)-YLim(1))*.285;
s = {['LSE: ' num2str(LSE)]};
text(xpos, ypos, s, 'color', 'r');
hold off
end



function f = plotErrBoundsWgh(X, XExt, Err, centiles, means, medians, mins, maxs, LSE, nSweeps, nameStr, type, expand)

f = figure('position', [50 50 1200 600]);
if strcmpi(type,'linear')
    p1 = plot(X, Err(1:length(X)), 'r');
    hold on
    p2 = plot(X, 0.9*Err(1:length(X)), 'r--');
    plot(X, 1.1*Err(1:length(X)), 'r--');
    plot(X, 0.8*Err(1:length(X)), 'r--');
    plot(X, 1.2*Err(1:length(X)), 'r--');
    plot(X, 0.7*Err(1:length(X)), 'r--');
    plot(X, 1.3*Err(1:length(X)), 'r--');
else
    p1 = plot(XExt, Err, 'r');
    hold on
    p2 = plot(XExt, [0.9*Err; 1.1*Err], 'r--');
    plot(XExt, [0.8*Err; 1.2*Err], 'r--');
    plot(XExt, [0.7*Err; 1.3*Err], 'r--');
end
p3 = plot(X, [mins'; maxs'], 'o', 'MarkerEdgeColor', 'y', 'MarkerFaceColor', 'y');
p4 = zeros(size(X));
p5 = p4;
p6 = p4;
for ix = 1:length(X)
    p4(ix) = plot(X(ix), medians(ix), 's', 'MarkerEdgeColor', 'b', 'MarkerFaceColor',...
        [1-((1-nSweeps(ix)/max(nSweeps))*1/max(nSweeps)) 1-nSweeps(ix)/max(nSweeps) 1-nSweeps(ix)/max(nSweeps)]);
    p5(ix) = plot(X(ix), means(ix), '^', 'MarkerEdgeColor', 'b', 'MarkerFaceColor',...
        [1-((1-nSweeps(ix)/max(nSweeps))*1/max(nSweeps)) 1-nSweeps(ix)/max(nSweeps) 1-nSweeps(ix)/max(nSweeps)]);
    p6(ix) = plot(X(ix), centiles(ix), 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor',...
        [1-((1-nSweeps(ix)/max(nSweeps))*1/max(nSweeps)) 1-nSweeps(ix)/max(nSweeps) 1-nSweeps(ix)/max(nSweeps)]);
end

set(f, 'Name', nameStr);
title(nameStr);
xlabel('Number of sweeps');
ylabel('Event count');
ind = (floor(ix/2));
if expand
    legend([p1 p2(1) p6(ind) p5(ind) p4(ind) p3(1)], 'Model fit to centiles', '10%,20%,30% deviation', '6-score combined 50th centile', 'Mean', 'Median',...
        'Range', 'Location', 'NorthWest');
else
    legend([p1 p2(1) p6(ind) p5(ind) p4(ind) p3(1)], 'Model fit to centiles', '10%,20%,30% deviation', '3-score combined 50th centile', 'Mean', 'Median',...
        'Range', 'Location', 'NorthWest');
end

XLim = xlim;
xpos = XLim(1) + (XLim(2)-XLim(1))*.01;
YLim = ylim;
ypos = YLim(2) - (YLim(2)-YLim(1))*.285;
s = {['LSE: ' num2str(LSE)]};
text(xpos, ypos, s, 'color', 'r');
hold off
end



function [Regr, LSE, coef] = lineMinisUncons(X, data, weights)

coef = 0:0.01:100;
offset = -200:500;
bestFits = zeros(length(offset), length(X));
leastErr = zeros(length(offset), 1);
bestCoefs = leastErr;
bestOffsets = leastErr;
for iC = 1:length(offset)
    regrArray = coef'*X + offset(iC);
    errArray = (regrArray - repmat(data, length(coef), 1)).^2;
    errArrayWgh = errArray.*repmat(weights, length(coef), 1);
    errSum = sum(errArrayWgh,2);
    [leastErr(iC), iLeastErr] = min(errSum);
    bestFits(iC,:) = regrArray(iLeastErr,:);
    bestCoefs(iC) = coef(iLeastErr);
    bestOffsets(iC) = offset(iC);
end

[~, iLSE] = min(leastErr);
coef = bestCoefs(iLSE);
Regr = coef*(1:6*length(X)) + bestOffsets(iLSE);
LSE = sum((Regr(1:length(data)) - data).^2);
end



function predictedUF = predictComparison(target, noise, pow, lengthRatio, XAmpsErrWgh_1, XAmpsDevErrWgh_1, XAmpsNegErrWgh_1, XAmpsNegDevErrWgh_1,...
    XRTsErrWgh_1, XRTsDevErrWgh_1, XRTsNegErrWgh_1, XRTsNegDevErrWgh_1, XTwoDsErrWgh_1, XTwoDsDevErrWgh_1, XTwoDsNegErrWgh_1, XTwoDsNegDevErrWgh_1)

diffAmps = .5*XAmpsErrWgh_1*(target^pow + lengthRatio*noise^pow);
diffAmpsDev = .5*XAmpsDevErrWgh_1*(target^pow + lengthRatio*noise^pow);
diffAmpsNeg = .5*XAmpsNegErrWgh_1*(target^pow + lengthRatio*noise^pow);
diffAmpsNegDev = .5*XAmpsNegDevErrWgh_1*(target^pow + lengthRatio*noise^pow);
diffRTs = .5*XRTsErrWgh_1*(target^pow + lengthRatio*noise^pow);
diffRTsDev = .5*XRTsDevErrWgh_1*(target^pow + lengthRatio*noise^pow);
diffRTsNeg = .5*XRTsNegErrWgh_1*(target^pow + lengthRatio*noise^pow);
diffRTsNegDev = .5*XRTsNegDevErrWgh_1*(target^pow + lengthRatio*noise^pow);
diffTwoDs = .5*XTwoDsErrWgh_1*(target^pow + lengthRatio*noise^pow);
diffTwoDsDev = .5*XTwoDsDevErrWgh_1*(target^pow + lengthRatio*noise^pow);
diffTwoDsNeg = .5*XTwoDsNegErrWgh_1*(target^pow + lengthRatio*noise^pow);
diffTwoDsNegDev = .5*XTwoDsNegDevErrWgh_1*(target^pow + lengthRatio*noise^pow);
predictedUF = [diffAmps; diffAmpsDev; diffAmpsNeg; diffAmpsNegDev; diffRTs; diffRTsDev; diffRTsNeg; diffRTsNegDev; diffTwoDs;...
    diffTwoDsDev; diffTwoDsNeg; diffTwoDsNegDev];
end



function predictedUF = predictComparisonAmps(target, noise, pow, lengthRatio, XAmpsErrWgh_1, XAmpsDevErrWgh_1, XAmpsNegErrWgh_1, XAmpsNegDevErrWgh_1)

diffAmps = .5*XAmpsErrWgh_1*(target^pow + lengthRatio*noise^pow);
diffAmpsDev = .5*XAmpsDevErrWgh_1*(target^pow + lengthRatio*noise^pow);
diffAmpsNeg = .5*XAmpsNegErrWgh_1*(target^pow + lengthRatio*noise^pow);
diffAmpsNegDev = .5*XAmpsNegDevErrWgh_1*(target^pow + lengthRatio*noise^pow);
predictedUF = [diffAmps; diffAmpsDev; diffAmpsNeg; diffAmpsNegDev];
end



function [f1, f2] = plotData(amplitudeArray, riseTimeArray, Amps, RTs, RTtype, fNames, flFactor, type)

f1 = figure('position', [50 50 1200 600]);
hold on
Amp99prc = zeros(1,size(Amps,1));
pAmp = Amp99prc;
for iData = 1:size(Amps,1)
    Ampcs = cumsum(Amps(iData,:))/sum(Amps(iData,:));
    Amp99prc(iData) = find(Ampcs > .99, 1);
end
iEndAmp = min([max(Amp99prc)+1 length(amplitudeArray)]);
xlim([amplitudeArray(1) amplitudeArray(iEndAmp)]);
colour = ([ones(size(Amps,1),1) zeros(size(Amps,1),2)]).*repmat((1/size(Amps,1): 1/size(Amps,1) :1)',1,3);
for iData = 1:size(Amps,1)
    pAmp(iData) = plot(amplitudeArray, flFactor(iData)*Amps(iData,:), '.-', 'color', colour(iData,:));
end
set(f1, 'NumberTitle', 'off');
set(f1, 'Name', 'One-dimensional Amplitude Histograms of Data Files');
title([type ' amplitude distributions']);
xlabel('Amplitude(mV)');
ylabel('Number of Events');
legend(pAmp, fNames, 'Location', 'NorthEast');
hold off

f2 = figure('position', [50 50 1200 600]);
hold on
RT99prc = zeros(1,size(RTs,1));
pRT = RT99prc;
for iData = 1:size(RTs,1)
    RTcs = cumsum(RTs(iData,:))/sum(RTs(iData,:));
    RT99prc(iData) = find(RTcs > .99, 1);
end
iEndRT = min([max(RT99prc)+1 length(riseTimeArray)]);
xlim([riseTimeArray(1) riseTimeArray(iEndRT)]);
for iData = 1:size(RTs,1)
    pRT(iData) = plot(riseTimeArray, flFactor(iData)*RTs(iData,:), '.-', 'color', colour(iData,:));
end
set(f2, 'NumberTitle', 'off');
set(f2, 'Name', 'One-dimensional Rise Time Histograms of Data Files');
if strcmpi(RTtype, '10-90%')
    title([type ' 10-90% rise time distributions']);
    xlabel('10-90% rise times (ms)');
elseif strcmpi(RTtype, '20-80%')
    title([type ' 20-80% rise time distributions']);
    xlabel('20-80% rise times (ms)');
end
ylabel('Number of Events');
legend(pRT, fNames, 'Location', 'NorthEast');
hold off
end



function f = plotCentiles(meanAmps, meanRTs, medianAmps, medianRTs, diffAmps, diffRTs, prctAmps, prctRTs, Name)

f = figure('position', [50 50 1200 600]);
dataToPlot = sortrows([diffAmps diffRTs]);
p1 = plot(dataToPlot(:,1)', dataToPlot(:,2)', '.r', 'MarkerSize', 5);
hold on
p2 = plot([min(dataToPlot(:,1)) meanAmps], [meanRTs meanRTs], '--k');
plot([meanAmps meanAmps], [min(dataToPlot(:,2)) meanRTs], '--k');
p3 = plot([min(dataToPlot(:,1)) medianAmps], [medianRTs medianRTs], '--g');
plot([medianAmps medianAmps], [min(dataToPlot(:,2)) medianRTs], '--g');
p4 = plot([min(dataToPlot(:,1)) prctAmps], [prctRTs prctRTs], '--b');
plot([prctAmps prctAmps], [min(dataToPlot(:,2)) prctRTs], '--b');
set(f, 'NumberTitle', 'off');
set(f, 'Name', 'Estimated 50th centile');
title(Name);
xlabel('Amplitude SAD');
ylabel('Rise time SAD');
legend([p1 p2 p3 p4],'Single sweep data', 'Mean', 'Median', 'Combined 50th centile', 'Location', 'NorthWest');
hold off
end



function prct = calc50thCentile(diffAmps, diffRTs, diffTwoDs, prct)

proceed = false;
while ~proceed
    diffAmpsPrct = diffAmps;
    diffRTsPrct = diffRTs;
    diffTwoDsPrct = diffTwoDs;
    
    prctAmps = prctile(diffAmpsPrct, prct);
    prctRTs = prctile(diffRTsPrct, prct);
    prctTwoDs = prctile(diffTwoDsPrct, prct);
    
    diffAmpsPrct(diffAmpsPrct > prctAmps) = 0;
    diffRTsPrct(diffRTsPrct > prctRTs) = 0;
    diffTwoDsPrct(diffTwoDsPrct > prctTwoDs) = 0;
    
    diffAmpsPrct = logical(diffAmpsPrct);
    diffRTsPrct = logical(diffRTsPrct);
    diffTwoDsPrct = logical(diffTwoDsPrct);
    
    combinedPrct = sum([diffAmpsPrct diffRTsPrct diffTwoDsPrct], 2);
    combinedPrct(combinedPrct < 3) = 0;
    combinedPrct = logical(combinedPrct);
    if sum(combinedPrct) < ceil(length(combinedPrct)/2)
        prct = prct + 1;
    else
        proceed = true;
    end
end
end



function prct = calc50thCentileExpand(diffAmps, diffAmpsBottom, diffAmpsMid, diffAmpsTop, diffRTs, diffTwoDs, prct)

proceed = false;
while ~proceed
    diffAmpsPrct = diffAmps;
    diffAmpsPrctBottom = diffAmpsBottom;
    diffAmpsPrctMid = diffAmpsMid;
    diffAmpsPrctTop = diffAmpsTop;
    diffRTsPrct = diffRTs;
    diffTwoDsPrct = diffTwoDs;
    
    prctAmps = prctile(diffAmpsPrct, prct);
    prctAmpsBottom = prctile(diffAmpsPrctBottom, prct);
    prctAmpsMid = prctile(diffAmpsPrctMid, prct);
    prctAmpsTop = prctile(diffAmpsPrctTop, prct);
    prctRTs = prctile(diffRTsPrct, prct);
    prctTwoDs = prctile(diffTwoDsPrct, prct);
    
    diffAmpsPrct(diffAmpsPrct > prctAmps) = 0;
    diffAmpsPrctBottom(diffAmpsPrctBottom > prctAmpsBottom) = 0;
    diffAmpsPrctMid(diffAmpsPrctMid > prctAmpsMid) = 0;
    diffAmpsPrctTop(diffAmpsPrctTop > prctAmpsTop) = 0;
    diffRTsPrct(diffRTsPrct > prctRTs) = 0;
    diffTwoDsPrct(diffTwoDsPrct > prctTwoDs) = 0;
    
    diffAmpsPrct = logical(diffAmpsPrct);
    diffAmpsPrctBottom = logical(diffAmpsPrctBottom);
    diffAmpsPrctMid = logical(diffAmpsPrctMid);
    diffAmpsPrctTop = logical(diffAmpsPrctTop);
    diffRTsPrct = logical(diffRTsPrct);
    diffTwoDsPrct = logical(diffTwoDsPrct);
    
    combinedPrct = sum([diffAmpsPrct diffAmpsPrctBottom diffAmpsPrctMid diffAmpsPrctTop diffRTsPrct diffTwoDsPrct], 2);
    combinedPrct(combinedPrct < 6) = 0;
    combinedPrct = logical(combinedPrct);
    if sum(combinedPrct) < ceil(length(combinedPrct)/2)
        prct = prct + 1;
    else
        proceed = true;
    end
end
end