function [diffAmpsAbs, diffAmpsNegAbs, diffAmps, diffAmpsNeg, diffAmpsMaxAbs, diffAmpsNegMaxAbs] = compareSweepsPartial(reduceComp, Amps, AmpsNeg)

compTable = nchoosek(1:size(Amps,1),2);
if reduceComp
    compTable(end,:) = [];
end
diffAmpsAbs = zeros(size(compTable,1),1);
diffAmpsNegAbs = diffAmpsAbs;
diffAmps = zeros(size(compTable,1),1);
diffAmpsNeg = diffAmps;
for n2s = 1:size(compTable,1)
    diffAmpsAbs(n2s) = sum(abs(Amps(compTable(n2s,1),:) - Amps(compTable(n2s,2),:)));
    diffAmpsNegAbs(n2s) = sum(abs(AmpsNeg(compTable(n2s,1),:) - AmpsNeg(compTable(n2s,2),:)));
    diffAmpsMaxAbs(n2s) = max(abs(Amps(compTable(n2s,1),:) - Amps(compTable(n2s,2),:)));
    diffAmpsNegMaxAbs(n2s) = max(abs(AmpsNeg(compTable(n2s,1),:) - AmpsNeg(compTable(n2s,2),:)));
    diffAmps(n2s) = sum(Amps(compTable(n2s,1),:) - Amps(compTable(n2s,2),:));
    diffAmpsNeg(n2s) = sum(AmpsNeg(compTable(n2s,1),:) - AmpsNeg(compTable(n2s,2),:));
end
end