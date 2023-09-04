function [noiseEst, dataProperties, F] = estimateNoiseMinis(filename, excludedTimes, targetMinis, detectionParameters, parallelCores, filtering)

F = zeros(1,2);
filenameSplit = regexp(filename, '\\*', 'split');
filenameShort = char(filenameSplit(end));

% Load an abf data file:
numSweeps = excludedTimes.nSweeps;
dataProperties = loadABF(filename);
if strcmpi(filtering.state,'on')
    [dataProperties.sweep, ~, F(1)] = filterMinis(dataProperties.sweep, dataProperties.dt, filtering, 'on');
end
detectionParameters.sampleInterval = dataProperties.dt;
detectionParameters.smoothWindow = round(detectionParameters.smoothWindow/detectionParameters.sampleInterval);
dataProperties.sweep = gsmooth(dataProperties.sweep, detectionParameters.smoothWindow, 3);

% Determine excluded times:
sweepDuration = (length(dataProperties.sweep)*dataProperties.dt - dataProperties.dt)/numSweeps;
startPulse = excludedTimes.startPulse;
endPulse = excludedTimes.endPulse;
startGlitch = excludedTimes.startGlitch;
endGlitch = excludedTimes.endGlitch;
excludedTimes = calcExcludedTimes(sweepDuration, numSweeps, 1000*startPulse, 1000*endPulse, 1000*startGlitch, 1000*endGlitch, dataProperties.dt);
excludedI = round(excludedTimes/dataProperties.dt) +1;

% Measure noise:
targetMinis(:,[1:2, 4:7, 10:end]) = [];
noiseEst = estimateNoise(dataProperties.sweep, excludedI, targetMinis, parallelCores);

% Mark noise estimates:
t = dataProperties.dt*(1:length(dataProperties.sweep)) - dataProperties.dt;
nameStringF = 'Summary: Noise estimates';
titleStringF = sprintf('Smoothed data from %s', filenameShort);
optionsF = struct('nameString', nameStringF, 'titleString', titleStringF, 'dataType', 'Membrane potential or current', 'dataUnits', '(mV or nA)');
F(2) = plotData(t, dataProperties.sweep, optionsF);
figure(F(2));
hold on;
plot(t(noiseEst(:,1)), noiseEst(:,4), '^', 'markersize', 10, 'markeredgecolor', 'r');
plot(t(noiseEst(:,2)), noiseEst(:,5), 'c.', 'markersize', 10);
plot(t(noiseEst(:,3)), noiseEst(:,5), 'b.', 'markersize', 10);
% Plot excluded times:
vexcluded = dataProperties.sweep(excludedI);
plot(excludedTimes, vexcluded, 'y.', 'markersize', 5);
hold off;
end


function noiseEst = estimateNoise(V, excludedI, targetMinis, parallelCores)

noiseMarkers = targetMinis;
row1 = zeros(1,length(excludedI));
for et = 1:length(excludedI)
    [row,~] = find(noiseMarkers == excludedI(et));
    if ~isempty(row)
        row1(et) = row;
    end
end
row1(row1==0) = [];
[row2,~] = find(noiseMarkers > length(V));
noiseMarkers(row2:end,:) = [];
noiseMarkers(row1,:) = [];
baselines = zeros(size(noiseMarkers,1),1);
tic
if parallelCores > 1
    parfor nm = 1:size(noiseMarkers,1)
        indices = noiseMarkers(nm,2:3); %#ok<*PFBNS>
        indices = indices(1):indices(2);
        baselines(nm) = mean(V(indices));
    end
else
    for nm = 1:size(noiseMarkers,1)
        baselines(nm) = mean(V(noiseMarkers(nm,2):noiseMarkers(nm,3)));
    end
end
toc
amplitudes = V(noiseMarkers(:,1))' - baselines;
noiseEst = [noiseMarkers V(noiseMarkers(:,1))' baselines amplitudes];
end