function [AmpsMean, AmpsMedian, AmpsMin, AmpsMax, AmpsPrct, meanAmpsDiff, meanAmpsDevDiff, meanAmpsNegDiff, meanAmpsNegDevDiff, medianAmpsDiff, medianAmpsDevDiff,...
    medianAmpsNegDiff, medianAmpsNegDevDiff, prctAmpsDiff, prctAmpsDevDiff, prctAmpsNegDiff, prctAmpsNegDevDiff, minAmpsDiff, minAmpsDevDiff, minAmpsNegDiff,...
    minAmpsNegDevDiff, maxAmpsDiff, maxAmpsDevDiff, maxAmpsNegDiff, maxAmpsNegDevDiff, initAmps, initAmpsNeg] = errorAmps(sweepcount, initAmps, initAmpsNeg,...
    prct, prctNeg, prctMAD, prctNegMAD, prctHisto)

diffAmps = struct([]);
diffAmpsDev = diffAmps;
diffAmpsNeg = diffAmps;
diffAmpsNegDev = diffAmps;
%
largestSweep = floor(sweepcount/2);
meanAmpsDiff = zeros(largestSweep,1);
meanAmpsDevDiff = meanAmpsDiff;
meanAmpsNegDiff = meanAmpsDiff;
meanAmpsNegDevDiff = meanAmpsDiff;
%
medianAmpsDiff = meanAmpsDiff;
medianAmpsDevDiff = meanAmpsDiff;
medianAmpsNegDiff = meanAmpsDiff;
medianAmpsNegDevDiff = meanAmpsDiff;
%
prctAmpsDiff = meanAmpsDiff;
prctAmpsDevDiff = meanAmpsDiff;
prctAmpsNegDiff = meanAmpsDiff;
prctAmpsNegDevDiff = meanAmpsDiff;
%
minAmpsDiff = meanAmpsDiff;
minAmpsDevDiff = meanAmpsDiff;
minAmpsNegDiff = meanAmpsDiff;
minAmpsNegDevDiff = meanAmpsDiff;
%
maxAmpsDiff = meanAmpsDiff;
maxAmpsDevDiff = meanAmpsDiff;
maxAmpsNegDiff = meanAmpsDiff;
maxAmpsNegDevDiff = meanAmpsDiff;
%
combAmps = sum(initAmps,1);
Ampcs = cumsum(combAmps)/sum(combAmps);
Ampcs = find(Ampcs > prctHisto, 1);
initAmps(:,1:Ampcs-1) = zeros(size(initAmps,1), Ampcs-1);
combAmpsNeg = sum(initAmpsNeg,1);
AmpcsNeg = cumsum(combAmpsNeg)/sum(combAmpsNeg);
AmpcsNeg = find(AmpcsNeg > prctHisto, 1);
initAmpsNeg(:,1:AmpcsNeg-1) = zeros(size(initAmpsNeg,1), AmpcsNeg-1);
%
nCompare = zeros(1,largestSweep);
nSweeps = ceil(sweepcount./(1:largestSweep));
reduceComp = nSweeps - floor(sweepcount./(1:largestSweep));
for sweepSize = 1:largestSweep
    nCompare(sweepSize) = nchoosek(nSweeps(sweepSize),2);
    if reduceComp(sweepSize)
        nCompare(sweepSize) = nCompare(sweepSize) - 1;
    end
    Amps = zeros(nSweeps(sweepSize),size(initAmps,2));
    AmpsNeg = Amps;
    for newSweep = 1:nSweeps(sweepSize)
        if newSweep ~= nSweeps(sweepSize)
            Amps(newSweep,:) = sum(initAmps((newSweep-1)*sweepSize+1: newSweep*sweepSize,:),1);
            AmpsNeg(newSweep,:) = sum(initAmpsNeg((newSweep-1)*sweepSize+1: newSweep*sweepSize,:),1);
        else
            Amps(newSweep,:) = sum(initAmps(end-sweepSize+1: end,:),1);
            AmpsNeg(newSweep,:) = sum(initAmpsNeg(end-sweepSize+1: end,:),1);
        end
    end
    [diffAmps{sweepSize}, diffAmpsDev{sweepSize}, diffAmpsNeg{sweepSize}, diffAmpsNegDev{sweepSize}] = compareSweepsAmps(reduceComp(sweepSize), Amps, AmpsNeg);
    
    meanAmpsDiff(sweepSize) = mean(diffAmps{sweepSize});
    meanAmpsDevDiff(sweepSize) = mean(diffAmpsDev{sweepSize});
    meanAmpsNegDiff(sweepSize) = mean(diffAmpsNeg{sweepSize});
    meanAmpsNegDevDiff(sweepSize) = mean(diffAmpsNegDev{sweepSize});
    
    medianAmpsDiff(sweepSize) = median(diffAmps{sweepSize});
    medianAmpsDevDiff(sweepSize) = median(diffAmpsDev{sweepSize});
    medianAmpsNegDiff(sweepSize) = median(diffAmpsNeg{sweepSize});
    medianAmpsNegDevDiff(sweepSize) = median(diffAmpsNegDev{sweepSize});
    
    minAmpsDiff(sweepSize) = min(diffAmps{sweepSize});
    minAmpsDevDiff(sweepSize) = min(diffAmpsDev{sweepSize});
    minAmpsNegDiff(sweepSize) = min(diffAmpsNeg{sweepSize});
    minAmpsNegDevDiff(sweepSize) = min(diffAmpsNegDev{sweepSize});
    
    maxAmpsDiff(sweepSize) = max(diffAmps{sweepSize});
    maxAmpsDevDiff(sweepSize) = max(diffAmpsDev{sweepSize});
    maxAmpsNegDiff(sweepSize) = max(diffAmpsNeg{sweepSize});
    maxAmpsNegDevDiff(sweepSize) = max(diffAmpsNegDev{sweepSize});
    
    prctAmpsDiff(sweepSize) = prctile(diffAmps{sweepSize},prct);
    prctAmpsDevDiff(sweepSize) = prctile(diffAmpsDev{sweepSize},prctMAD);
    prctAmpsNegDiff(sweepSize) = prctile(diffAmpsNeg{sweepSize},prctNeg);
    prctAmpsNegDevDiff(sweepSize) = prctile(diffAmpsNegDev{sweepSize},prctNegMAD);
end
X = (1:sweepSize);
AmpsMean = [X; meanAmpsDiff'; meanAmpsDevDiff'; meanAmpsNegDiff'; meanAmpsNegDevDiff'];
AmpsMedian = [X; medianAmpsDiff'; medianAmpsDevDiff'; medianAmpsNegDiff'; medianAmpsNegDevDiff'];
AmpsMin = [X; minAmpsDiff'; minAmpsDevDiff'; minAmpsNegDiff'; minAmpsNegDevDiff'];
AmpsMax = [X; maxAmpsDiff'; maxAmpsDevDiff'; maxAmpsNegDiff'; maxAmpsNegDevDiff'];
AmpsPrct = [X; prctAmpsDiff'; prctAmpsDevDiff'; prctAmpsNegDiff'; prctAmpsNegDevDiff'];
end





%% Local functions:
function [diffAmps, diffAmpsDev, diffAmpsNeg, diffAmpsNegDev] = compareSweepsAmps(reduceComp, Amps, AmpsNeg)

compTable = nchoosek(1:size(Amps,1),2);
if reduceComp
    compTable(end,:) = [];
end
diffAmps = zeros(size(compTable,1),1);
diffAmpsNeg = diffAmps;
diffAmpsDev = diffAmps;
diffAmpsNegDev = diffAmps;
for n2s = 1:size(compTable,1)
    diffAmps(n2s) = sum(abs(Amps(compTable(n2s,1),:) - Amps(compTable(n2s,2),:)));
    diffAmpsDev(n2s) = max(abs(Amps(compTable(n2s,1),:) - Amps(compTable(n2s,2),:)));
    diffAmpsNeg(n2s) = sum(abs(AmpsNeg(compTable(n2s,1),:) - AmpsNeg(compTable(n2s,2),:)));
    diffAmpsNegDev(n2s) = max(abs(AmpsNeg(compTable(n2s,1),:) - AmpsNeg(compTable(n2s,2),:)));
end
end