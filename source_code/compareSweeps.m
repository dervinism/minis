function [diffAmps, diffAmpsDev, diffAmpsNeg, diffAmpsNegDev, diffRTs, diffRTsDev, diffRTsNeg, diffRTsNegDev, diff2Ds, diff2DsDev, diff2DsNeg,...
    diff2DsNegDev] = compareSweeps(reduceComp, Amps, AmpsNeg, RTs, RTsNeg, twoDs, twoDsNeg)

compTable = nchoosek(1:size(Amps,1),2);
if reduceComp
    compTable(end,:) = [];
end
diffAmps = zeros(size(compTable,1),1);
diffAmpsNeg = diffAmps;
diffRTs = diffAmps;
diffRTsNeg = diffAmps;
diff2Ds = diffAmps;
diff2DsNeg = diffAmps;
diffAmpsDev = diffAmps;
diffAmpsNegDev = diffAmps;
diffRTsDev = diffAmps;
diffRTsNegDev = diffAmps;
diff2DsDev = diffAmps;
diff2DsNegDev = diffAmps;
for n2s = 1:size(compTable,1)
    diffAmps(n2s) = sum(abs(Amps(compTable(n2s,1),:) - Amps(compTable(n2s,2),:)));
    diffAmpsDev(n2s) = max(abs(Amps(compTable(n2s,1),:) - Amps(compTable(n2s,2),:)));
    diffAmpsNeg(n2s) = sum(abs(AmpsNeg(compTable(n2s,1),:) - AmpsNeg(compTable(n2s,2),:)));
    diffAmpsNegDev(n2s) = max(abs(AmpsNeg(compTable(n2s,1),:) - AmpsNeg(compTable(n2s,2),:)));
    diffRTs(n2s) = sum(abs(RTs(compTable(n2s,1),:) - RTs(compTable(n2s,2),:)));
    diffRTsDev(n2s) = max(abs(RTs(compTable(n2s,1),:) - RTs(compTable(n2s,2),:)));
    diffRTsNeg(n2s) = sum(abs(RTsNeg(compTable(n2s,1),:) - RTsNeg(compTable(n2s,2),:)));
    diffRTsNegDev(n2s) = max(abs(RTsNeg(compTable(n2s,1),:) - RTsNeg(compTable(n2s,2),:)));
    diff2Ds(n2s) = sum(sum(abs(twoDs{compTable(n2s,1)} - twoDs{compTable(n2s,2)})));
    diff2DsDev(n2s) = max(max(abs(twoDs{compTable(n2s,1)} - twoDs{compTable(n2s,2)})));
    diff2DsNeg(n2s) = sum(sum(abs(twoDsNeg{compTable(n2s,1)} - twoDsNeg{compTable(n2s,2)})));
    diff2DsNegDev(n2s) = max(max(abs(twoDsNeg{compTable(n2s,1)} - twoDsNeg{compTable(n2s,2)})));
end
end