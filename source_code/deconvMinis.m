function F = deconvMinis(targetFilename,noiseFilename,targetExcludedTimes,noiseExcludedTimes,detectionParameters,classificationParameters,filtering,parallelCores)

F = zeros(1,6);
if ~isempty(noiseExcludedTimes.endPulse)
    noiseFilt = struct('state', filtering.state, 'nSweeps', filtering.nNoiseSweeps, 'pulseEnd', noiseExcludedTimes.endPulse(end));
else
    noiseFilt = struct('state', filtering.state, 'nSweeps', filtering.nNoiseSweeps, 'pulseEnd', 0);
end
if ~isempty(targetExcludedTimes.endPulse)
    targetFilt = struct('state', filtering.state, 'nSweeps', filtering.nTargetSweeps, 'pulseEnd', targetExcludedTimes.endPulse(end));
else
    targetFilt = struct('state', filtering.state, 'nSweeps', filtering.nTargetSweeps, 'pulseEnd', 0);
end

[targetMinis, ~, ~, targetProperties, F(1:2)] = detectMinisStandalone('off', targetFilename, detectionParameters, targetExcludedTimes,...
    classificationParameters, parallelCores, 'none', 1, 0, noiseFilt);
[noiseEstimates, noiseProperties, F(3:4)] = estimateNoiseMinis(noiseFilename,noiseExcludedTimes,targetMinis,detectionParameters,parallelCores,targetFilt);

amplitudeArray = classificationParameters.amplitudeArray;
amplitudeArrayExt = classificationParameters.amplitudeArrayExt;

lengthRatio = (length(targetProperties.sweep)-length(targetExcludedTimes))/(length(noiseProperties.sweep)-length(noiseExcludedTimes));
targetAmp = hist(targetMinis(:,4), amplitudeArrayExt); %#ok<*HIST>
targetAmp = targetAmp(1:end-1);
ampstep = amplitudeArray(2) - amplitudeArray(1);
if amplitudeArrayExt(1) == 0
    noiseAmpArrayExt = [-fliplr(amplitudeArrayExt) amplitudeArrayExt(2:end)];
else
    noiseAmpArrayExt = [-fliplr(amplitudeArrayExt) (-amplitudeArrayExt(1)+ampstep: ampstep :amplitudeArray(1)-ampstep) amplitudeArrayExt];
end
noiseAmp = hist(noiseEstimates(:,6), noiseAmpArrayExt);
noiseAmpArray = noiseAmpArrayExt(2:end-1);
noiseAmp = lengthRatio*noiseAmp(2:end-1);
targetAmp = [zeros(1, length(noiseAmp) - length(targetAmp)) targetAmp];
[~, minis] = deconv(targetAmp+1, noiseAmp+1);
% Fna = fft([noiseAmp(ceil(length(noiseAmp)/2):end) noiseAmp(1:floor(length(noiseAmp)/2))]);
% Fta = fft([targetAmp(ceil(length(targetAmp)/2):end) targetAmp(1:floor(length(targetAmp)/2))]);
% Fminis = Fna/Fta;
% minis = ifft(Fminis);
F(5) = plotConv(amplitudeArray, noiseAmpArray, targetAmp, noiseAmp, minis);

function f = plotConv(amplitudeArray, noiseAmpArray, targetAmp, noiseAmp, minis)
f = figure('position', [50 50 1200 600]);
figure(f);
hold on
iEndAmp = min([find(targetAmp,1,'last')+1 length(amplitudeArray)]);
xlim([-amplitudeArray(iEndAmp) amplitudeArray(iEndAmp)]);
plot(noiseAmpArray, targetAmp, 'b.-');
plot(noiseAmpArray, minis, 'g.-');
plot(noiseAmpArray, noiseAmp, 'k.-');
set(f, 'NumberTitle', 'off');
set(f, 'Name', 'One-dimensional histograms');
title('Amplitude distribution');
xlabel('Amplitude(mV)');
ylabel('Number of Events');
hold off