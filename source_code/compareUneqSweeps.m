function [estimated, prctAmps, prctAmpsDev, prctAmpsNeg, prctAmpsNegDev, prctRTs, prctRTsDev, prctRTsNeg, prctRTsNegDev, prctTwoDs, prctTwoDsDev, prctTwoDsNeg,...
    prctTwoDsNegDev, minAmps, minAmpsDev, minAmpsNeg, minAmpsNegDev, minRTs, minRTsDev, minRTsNeg, minRTsNegDev, minTwoDs, minTwoDsDev, minTwoDsNeg, minTwoDsNegDev,...
    maxAmps, maxAmpsDev, maxAmpsNeg, maxAmpsNegDev, maxRTs, maxRTsDev, maxRTsNeg, maxRTsNegDev, maxTwoDs, maxTwoDsDev, maxTwoDsNeg, maxTwoDsNegDev,...
    lengthRatio] = compareUneqSweeps(sweepcount, filtfs, Amps, AmpsNeg, RTs, RTsNeg, TwoDs, TwoDsNeg, prct, prctNeg, prctMAD, prctNegMAD)

target = str2double(filtfs{1});
noise = str2double(filtfs{2});
lengthRatio = target/noise;
if (target + noise) > sweepcount || target == noise || target < 1 || noise < 1
    estimated = false;
    prctAmps = [];
    prctAmpsDev = [];
    prctAmpsNeg = [];
    prctAmpsNegDev = [];
    prctRTs = [];
    prctRTsDev = [];
    prctRTsNeg = [];
    prctRTsNegDev = [];
    prctTwoDs = [];
    prctTwoDsDev = [];
    prctTwoDsNeg = [];
    prctTwoDsNegDev = [];
    minAmps = [];
    minAmpsDev = [];
    minAmpsNeg = [];
    minAmpsNegDev = [];
    minRTs = [];
    minRTsDev = [];
    minRTsNeg = [];
    minRTsNegDev = [];
    minTwoDs = [];
    minTwoDsDev = [];
    minTwoDsNeg = [];
    minTwoDsNegDev = [];
    maxAmps = [];
    maxAmpsDev = [];
    maxAmpsNeg = [];
    maxAmpsNegDev = [];
    maxRTs = [];
    maxRTsDev = [];
    maxRTsNeg = [];
    maxRTsNegDev = [];
    maxTwoDs = [];
    maxTwoDsDev = [];
    maxTwoDsNeg = [];
    maxTwoDsNegDev = [];
    if (target + noise) > sweepcount
        uiwait(msgbox(...
            'Because the error bound cannot be estimated based on the data for the required file length, only a theoretical prediction will be provided.',...
            'modal'));
        return
    end
end

diffAmps = [];
diffAmpsNeg = [];
diffRTs = [];
diffRTsNeg = [];
diffTwoDs = [];
diffTwoDsNeg = [];
diffAmpsDev = [];
diffAmpsNegDev = [];
diffRTsDev = [];
diffRTsNegDev = [];
diffTwoDsDev = [];
diffTwoDsNegDev = [];

nTarget = ceil(sweepcount/target);
targEnd = sweepcount - (nTarget-1:-1:1)*target;
iTargEnd = find(targEnd >= noise, 1);
nTargRed = nTarget - iTargEnd;
for iTarget = 1:nTargRed
    if iTarget ~= nTarget
        AmpsTarg = Amps((iTarget-1)*target+1 : iTarget*target,:);
        AmpsNegTarg = AmpsNeg((iTarget-1)*target+1 : iTarget*target,:);
        RTsTarg = RTs((iTarget-1)*target+1 : iTarget*target,:);
        RTsNegTarg = RTsNeg((iTarget-1)*target+1 : iTarget*target,:);
        TwoDsTarg = TwoDs((iTarget-1)*target+1 : iTarget*target);
        TwoDsNegTarg = TwoDsNeg((iTarget-1)*target+1 : iTarget*target);
        AmpsNoise = Amps(iTarget*target+1:end,:);
        AmpsNegNoise = AmpsNeg(iTarget*target+1:end,:);
        RTsNoise = RTs(iTarget*target+1:end,:);
        RTsNegNoise = RTsNeg(iTarget*target+1:end,:);
        TwoDsNoise = TwoDs(iTarget*target+1:end);
        TwoDsNegNoise = TwoDsNeg(iTarget*target+1:end);
    else
        AmpsTarg = Amps(end-target+1 : end,:);
        AmpsNegTarg = AmpsNeg(end-target+1 : end,:);
        RTsTarg = RTs(end-target+1 : end,:);
        RTsNegTarg = RTsNeg(end-target+1 : end,:);
        TwoDsTarg = TwoDs(end-target+1 : end);
        TwoDsNegTarg = TwoDsNeg(end-target+1 : end);
        AmpsNoise = Amps(iTarget*target+1:end,:);
        AmpsNegNoise = AmpsNeg(iTarget*target+1:end,:);
        RTsNoise = RTs(iTarget*target+1:end,:);
        RTsNegNoise = RTsNeg(iTarget*target+1:end,:);
        TwoDsNoise = TwoDs(iTarget*target+1:end);
        TwoDsNegNoise = TwoDsNeg(iTarget*target+1:end);
    end
    
    nNoise = ceil(size(AmpsNoise,1)/noise);
    diffAmpsTemp = zeros(nNoise,1);
    diffAmpsDevTemp = diffAmpsTemp;
    diffAmpsNegTemp = diffAmpsTemp;
    diffAmpsNegDevTemp  = diffAmpsTemp;
    diffRTsTemp = diffAmpsTemp;
    diffRTsDevTemp = diffAmpsTemp;
    diffRTsNegTemp = diffAmpsTemp;
    diffRTsNegDevTemp = diffAmpsTemp;
    diffTwoDsTemp = diffAmpsTemp;
    diffTwoDsNegTemp = diffAmpsTemp;
    diffTwoDsDevTemp = diffAmpsTemp;
    diffTwoDsNegDevTemp = diffAmpsTemp;
    for iNoise = 1:nNoise
        if iNoise ~= nNoise
            diffAmpsTemp(iNoise) = sum(abs(sum(AmpsTarg,1) - lengthRatio*sum(AmpsNoise((iNoise-1)*noise+1:iNoise*noise,:),1)));
            diffAmpsDevTemp(iNoise) = max(abs(sum(AmpsTarg,1) - lengthRatio*sum(AmpsNoise((iNoise-1)*noise+1:iNoise*noise,:),1)));
            diffAmpsNegTemp(iNoise) = sum(abs(sum(AmpsNegTarg,1) - lengthRatio*sum(AmpsNegNoise((iNoise-1)*noise+1:iNoise*noise,:),1)));
            diffAmpsNegDevTemp(iNoise) = max(abs(sum(AmpsNegTarg,1) - lengthRatio*sum(AmpsNegNoise((iNoise-1)*noise+1:iNoise*noise,:),1)));
            diffRTsTemp(iNoise) = sum(abs(sum(RTsTarg,1) - lengthRatio*sum(RTsNoise((iNoise-1)*noise+1:iNoise*noise,:),1)));
            diffRTsDevTemp(iNoise) = max(abs(sum(RTsTarg,1) - lengthRatio*sum(RTsNoise((iNoise-1)*noise+1:iNoise*noise,:),1)));
            diffRTsNegTemp(iNoise) = sum(abs(sum(RTsNegTarg,1) - lengthRatio*sum(RTsNegNoise((iNoise-1)*noise+1:iNoise*noise,:),1)));
            diffRTsNegDevTemp(iNoise) = max(abs(sum(RTsNegTarg,1) - lengthRatio*sum(RTsNegNoise((iNoise-1)*noise+1:iNoise*noise,:),1)));
            targTemp = cell2mat(TwoDsTarg(1));
            targNegTemp = cell2mat(TwoDsNegTarg(1));
            for iCellT = 2:target
                targTemp = targTemp + cell2mat(TwoDsTarg(iCellT));
                targNegTemp = targNegTemp + cell2mat(TwoDsNegTarg(iCellT));
            end
            TwoDsNoiseTemp = cell2mat(TwoDsNoise((iNoise-1)*noise+1));
            TwoDsNegNoiseTemp = cell2mat(TwoDsNegNoise((iNoise-1)*noise+1));
            for iCellN = 2:noise
                TwoDsNoiseTemp = TwoDsNoiseTemp + cell2mat(TwoDsNoise((iNoise-1)*noise+iCellN));
                TwoDsNegNoiseTemp = TwoDsNegNoiseTemp + cell2mat(TwoDsNegNoise((iNoise-1)*noise+iCellN));
            end
            diffTwoDsTemp(iNoise) = sum(sum(abs(targTemp - lengthRatio*TwoDsNoiseTemp)));
            diffTwoDsDevTemp(iNoise) = max(max(abs(targTemp - lengthRatio*TwoDsNoiseTemp)));
            diffTwoDsNegTemp(iNoise) = sum(sum(abs(targNegTemp - lengthRatio*TwoDsNegNoiseTemp)));
            diffTwoDsNegDevTemp(iNoise) = max(max(abs(targNegTemp - lengthRatio*TwoDsNegNoiseTemp)));
        else
            diffAmpsTemp(iNoise) = sum(abs(sum(AmpsTarg,1) - lengthRatio*sum(AmpsNoise(end-noise+1:end,:),1)));
            diffAmpsDevTemp(iNoise) = max(abs(sum(AmpsTarg,1) - lengthRatio*sum(AmpsNoise(end-noise+1:end,:),1)));
            diffAmpsNegTemp(iNoise) = sum(abs(sum(AmpsNegTarg,1) - lengthRatio*sum(AmpsNegNoise(end-noise+1:end,:),1)));
            diffAmpsNegDevTemp(iNoise) = max(abs(sum(AmpsNegTarg,1) - lengthRatio*sum(AmpsNegNoise(end-noise+1:end,:),1)));
            diffRTsTemp(iNoise) = sum(abs(sum(RTsTarg,1) - lengthRatio*sum(RTsNoise(end-noise+1:end,:),1)));
            diffRTsDevTemp(iNoise) = max(abs(sum(RTsTarg,1) - lengthRatio*sum(RTsNoise(end-noise+1:end,:),1)));
            diffRTsNegTemp(iNoise) = sum(abs(sum(RTsNegTarg,1) - lengthRatio*sum(RTsNegNoise(end-noise+1:end,:),1)));
            diffRTsNegDevTemp(iNoise) = max(abs(sum(RTsNegTarg,1) - lengthRatio*sum(RTsNegNoise(end-noise+1:end,:),1)));
            targTemp = cell2mat(TwoDsTarg(1));
            targNegTemp = cell2mat(TwoDsNegTarg(1));
            for iCellT = 2:target
                targTemp = targTemp + cell2mat(TwoDsTarg(iCellT));
                targNegTemp = targNegTemp + cell2mat(TwoDsNegTarg(iCellT));
            end
            TwoDsNoiseTemp = cell2mat(TwoDsNoise(end-noise+1));
            TwoDsNegNoiseTemp = cell2mat(TwoDsNegNoise(end-noise+1));
            for iCellN = 2:noise
                TwoDsNoiseTemp = TwoDsNoiseTemp + cell2mat(TwoDsNoise(end-noise+iCellN));
                TwoDsNegNoiseTemp = TwoDsNegNoiseTemp + cell2mat(TwoDsNegNoise(end-noise+iCellN));
            end
            diffTwoDsTemp(iNoise) = sum(sum(abs(targTemp - lengthRatio*TwoDsNoiseTemp)));
            diffTwoDsDevTemp(iNoise) = max(max(abs(targTemp - lengthRatio*TwoDsNoiseTemp)));
            diffTwoDsNegTemp(iNoise) = sum(sum(abs(targNegTemp - lengthRatio*TwoDsNegNoiseTemp)));
            diffTwoDsNegDevTemp(iNoise) = max(max(abs(targNegTemp - lengthRatio*TwoDsNegNoiseTemp)));
        end
    end
    
    diffAmps = [diffAmps; diffAmpsTemp]; %#ok<*AGROW>
    diffAmpsDev = [diffAmpsDev; diffAmpsDevTemp];
    diffAmpsNeg = [diffAmpsNeg; diffAmpsNegTemp];
    diffAmpsNegDev = [diffAmpsNegDev; diffAmpsNegDevTemp];
    diffRTs = [diffRTs; diffRTsTemp];
    diffRTsDev = [diffRTsDev; diffRTsDevTemp];
    diffRTsNeg = [diffRTsNeg; diffRTsNegTemp];
    diffRTsNegDev = [diffRTsNegDev; diffRTsNegDevTemp];
    diffTwoDs = [diffTwoDs; diffTwoDsTemp];
    diffTwoDsDev = [diffTwoDsDev; diffTwoDsDevTemp];
    diffTwoDsNeg = [diffTwoDsNeg; diffTwoDsNegTemp];
    diffTwoDsNegDev = [diffTwoDsNegDev; diffTwoDsNegDevTemp];
end

prctAmps = prctile(diffAmps, prct);
prctAmpsDev = prctile(diffAmpsDev, prctMAD);
prctAmpsNeg = prctile(diffAmpsNeg, prctNeg);
prctAmpsNegDev = prctile(diffAmpsNegDev, prctNegMAD);
prctRTs = prctile(diffRTs, prct);
prctRTsDev = prctile(diffRTsDev, prctMAD);
prctRTsNeg = prctile(diffRTsNeg, prctNeg);
prctRTsNegDev = prctile(diffRTsNegDev, prctNegMAD);
prctTwoDs = prctile(diffTwoDs, prct);
prctTwoDsDev = prctile(diffTwoDsDev, prctMAD);
prctTwoDsNeg = prctile(diffTwoDsNeg, prctNeg);
prctTwoDsNegDev = prctile(diffTwoDsNegDev, prctNegMAD);

minAmps = min(diffAmps);
minAmpsDev = min(diffAmpsDev);
minAmpsNeg = min(diffAmpsNeg);
minAmpsNegDev = min(diffAmpsNegDev);
minRTs = min(diffRTs);
minRTsDev = min(diffRTsDev);
minRTsNeg = min(diffRTsNeg);
minRTsNegDev = min(diffRTsNegDev);
minTwoDs = min(diffTwoDs);
minTwoDsDev = min(diffTwoDsDev);
minTwoDsNeg = min(diffTwoDsNeg);
minTwoDsNegDev = min(diffTwoDsNegDev);

maxAmps = max(diffAmps);
maxAmpsDev = max(diffAmpsDev);
maxAmpsNeg = max(diffAmpsNeg);
maxAmpsNegDev = max(diffAmpsNegDev);
maxRTs = max(diffRTs);
maxRTsDev = max(diffRTsDev);
maxRTsNeg = max(diffRTsNeg);
maxRTsNegDev = max(diffRTsNegDev);
maxTwoDs = max(diffTwoDs);
maxTwoDsDev = max(diffTwoDsDev);
maxTwoDsNeg = max(diffTwoDsNeg);
maxTwoDsNegDev = max(diffTwoDsNegDev);

estimated = true;