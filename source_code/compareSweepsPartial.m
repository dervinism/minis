function [diffAmps, diffAmpsNeg] = compareSweepsPartial(reduceComp, Amps, AmpsNeg)

compTable = nchoosek(1:size(Amps,1),2);
if reduceComp
    compTable(end,:) = [];
end
diffAmps = zeros(size(compTable,1),1);
diffAmpsNeg = diffAmps;
for n2s = 1:size(compTable,1)
    diffAmps(n2s) = sum(abs(Amps(compTable(n2s,1),:) - Amps(compTable(n2s,2),:)));
    diffAmpsNeg(n2s) = sum(abs(AmpsNeg(compTable(n2s,1),:) - AmpsNeg(compTable(n2s,2),:)));
end
end