function [prctAmps, prctAmpsDev, prctAmpsNeg, prctAmpsNegDev, minAmps, minAmpsDev, minAmpsNeg, minAmpsNegDev, maxAmps, maxAmpsDev, maxAmpsNeg,...
    maxAmpsNegDev] = compareSweepsAmps(sweepcount, filtfs, Amps, AmpsNeg, prct, prctNeg, prctMAD, prctNegMAD)

target = str2double(filtfs{1});
noise = str2double(filtfs{2});
lengthRatio = target/noise;

diffAmps = [];
diffAmpsNeg = [];
diffAmpsDev = [];
diffAmpsNegDev = [];

nTarget = ceil(sweepcount/target);
targEnd = sweepcount - (nTarget-1:-1:1)*target;
iTargEnd = find(targEnd >= noise, 1);
nTargRed = nTarget - iTargEnd;
for iTarget = 1:nTargRed
    if iTarget ~= nTarget
        AmpsTarg = Amps((iTarget-1)*target+1 : iTarget*target,:);
        AmpsNegTarg = AmpsNeg((iTarget-1)*target+1 : iTarget*target,:);
        AmpsNoise = Amps(iTarget*target+1:end,:);
        AmpsNegNoise = AmpsNeg(iTarget*target+1:end,:);
    else
        AmpsTarg = Amps(end-target+1 : end,:);
        AmpsNegTarg = AmpsNeg(end-target+1 : end,:);
        AmpsNoise = Amps(iTarget*target+1:end,:);
        AmpsNegNoise = AmpsNeg(iTarget*target+1:end,:);
    end
    
    nNoise = ceil(size(AmpsNoise,1)/noise);
    diffAmpsTemp = zeros(nNoise,1);
    diffAmpsDevTemp = diffAmpsTemp;
    diffAmpsNegTemp = diffAmpsTemp;
    diffAmpsNegDevTemp  = diffAmpsTemp;
    for iNoise = 1:nNoise
        if iNoise ~= nNoise
            diffAmpsTemp(iNoise) = sum(abs(sum(AmpsTarg,1) - lengthRatio*sum(AmpsNoise((iNoise-1)*noise+1:iNoise*noise,:),1)));
            diffAmpsDevTemp(iNoise) = max(abs(sum(AmpsTarg,1) - lengthRatio*sum(AmpsNoise((iNoise-1)*noise+1:iNoise*noise,:),1)));
            diffAmpsNegTemp(iNoise) = sum(abs(sum(AmpsNegTarg,1) - lengthRatio*sum(AmpsNegNoise((iNoise-1)*noise+1:iNoise*noise,:),1)));
            diffAmpsNegDevTemp(iNoise) = max(abs(sum(AmpsNegTarg,1) - lengthRatio*sum(AmpsNegNoise((iNoise-1)*noise+1:iNoise*noise,:),1)));
        else
            diffAmpsTemp(iNoise) = sum(abs(sum(AmpsTarg,1) - lengthRatio*sum(AmpsNoise(end-noise+1:end,:),1)));
            diffAmpsDevTemp(iNoise) = max(abs(sum(AmpsTarg,1) - lengthRatio*sum(AmpsNoise(end-noise+1:end,:),1)));
            diffAmpsNegTemp(iNoise) = sum(abs(sum(AmpsNegTarg,1) - lengthRatio*sum(AmpsNegNoise(end-noise+1:end,:),1)));
            diffAmpsNegDevTemp(iNoise) = max(abs(sum(AmpsNegTarg,1) - lengthRatio*sum(AmpsNegNoise(end-noise+1:end,:),1)));
        end
    end
    
    diffAmps = [diffAmps; diffAmpsTemp]; %#ok<*AGROW>
    diffAmpsDev = [diffAmpsDev; diffAmpsDevTemp];
    diffAmpsNeg = [diffAmpsNeg; diffAmpsNegTemp];
    diffAmpsNegDev = [diffAmpsNegDev; diffAmpsNegDevTemp];
end

prctAmps = prctile(diffAmps, prct);
prctAmpsDev = prctile(diffAmpsDev, prctMAD);
prctAmpsNeg = prctile(diffAmpsNeg, prctNeg);
prctAmpsNegDev = prctile(diffAmpsNegDev, prctNegMAD);

minAmps = min(diffAmps);
minAmpsDev = min(diffAmpsDev);
minAmpsNeg = min(diffAmpsNeg);
minAmpsNegDev = min(diffAmpsNegDev);

maxAmps = max(diffAmps);
maxAmpsDev = max(diffAmpsDev);
maxAmpsNeg = max(diffAmpsNeg);
maxAmpsNegDev = max(diffAmpsNegDev);