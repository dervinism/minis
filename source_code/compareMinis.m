function [targetMinis, noiseMinis, targetWave, F] = compareMinis(targetFilename, noiseFilename, targetExcludedTimes, noiseExcludedTimes, detectionParameters,...
    classificationParameters, filtering, waveform, parallelCores)
% COMPAREMINIS is helper function of MINIS GUI. For more information see
% the code.
%

F = zeros(1,8);
waveformT.estimate = 1;
waveformT.tau_m = waveform.tau_m;
waveformT.riseTimeArrayExt = waveform.riseTimeArrayExt;
waveformN.estimate = false;
noiseFilt = struct('state', filtering.state, 'excludedTimes', noiseExcludedTimes);
targetFilt = struct('state', filtering.state, 'excludedTimes', targetExcludedTimes);

[targetMinis, ~, targetWave, targetProperties, F(1:2), targetSpectra] = detectMinisStandalone(false, targetFilename, detectionParameters,...
    targetExcludedTimes, classificationParameters, parallelCores, 'none', true, waveformT, targetFilt, [1 1 1 1]);
close(F(2));
if isempty(targetSpectra)
    targetFilt.state = 'spectrum';
    targetFilt.nSweeps = targetProperties.hd.lActualEpisodes;
    [~, targetSpectra] = filterMinis(targetProperties.sweep, targetProperties.dt, targetFilt, false);
end
if isfield(targetWave,'parameters')
    targetMinis = [targetMinis repmat([targetWave.parameters.averageAmp targetWave.parameters.medianAmp targetWave.parameters.tau_m], size(targetMinis,1),1)];
end
[noiseMinis, ~, ~, noiseProperties, F(3:4), noiseSpectra] = detectMinisStandalone(false, noiseFilename, detectionParameters, noiseExcludedTimes,...
    classificationParameters, parallelCores, 'none', true, waveformN, noiseFilt, [1 1 1 0]);
close(F(4));
if isempty(noiseSpectra)
    noiseFilt.state = 'spectrum';
    noiseFilt.nSweeps = noiseProperties.hd.lActualEpisodes;
    [~, noiseSpectra] = filterMinis(noiseProperties.sweep, noiseProperties.dt, noiseFilt, false);
end

lengthRatio = (length(targetProperties.sweep)-length(targetExcludedTimes))/(length(noiseProperties.sweep)-length(noiseExcludedTimes));
[~, ~, target2D] = classifyMinis(targetMinis(:,4), targetMinis(:,12), classificationParameters);
[~, ~, noise2D] = classifyMinis(noiseMinis(:,4), noiseMinis(:,12), classificationParameters);
noise2D = lengthRatio*noise2D;
[F(5) F(6) F(7) F(8)] = plotCompare(classificationParameters.amplitudeArray, classificationParameters.riseTimeArray, target2D, noise2D,...
    detectionParameters.RTinterval, targetSpectra(2:3,:), noiseSpectra(2:3,:), noiseSpectra(4,:));