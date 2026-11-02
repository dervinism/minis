function [targetMinis, noiseMinis, targetWave, F, targetIncidenceRates, noiseIncidenceRates] = compareMinis(targetFilename, noiseFilename, targetExcludedTimes, noiseExcludedTimes, detectionParameters,...
    classificationParameters, filtering, waveform, parallelCores, showProgress)
% COMPAREMINIS is helper function of MINIS GUI. For more information see
% the code.
%
% SHOWPROGRESS (optional, default false) displays a progress bar during
% each of the two detection passes (target then noise). It is opt-in and
% should only be enabled by the GUI's Detect-and-compare task.

if nargin < 10 || isempty(showProgress)
    showProgress = false;
end
if showProgress
    targetProgress = struct('enable', true, 'label', 'Detecting in target file');
    noiseProgress = struct('enable', true, 'label', 'Detecting in noise file');
else
    targetProgress = struct('enable', false, 'label', '');
    noiseProgress = struct('enable', false, 'label', '');
end

F = zeros(1,8);
waveformT.estimate = 1;
waveformT.tau_m = waveform.tau_m;
waveformT.riseTimeArrayExt = waveform.riseTimeArrayExt;
waveformN.estimate = false;
noiseFilt = struct('state', filtering.state, 'excludedTimes', noiseExcludedTimes);
targetFilt = struct('state', filtering.state, 'excludedTimes', targetExcludedTimes);

[targetMinis, ~, targetWave, targetProperties, F(1:2), targetSpectra] = detectMinisStandalone(false, targetFilename, detectionParameters,...
    targetExcludedTimes, classificationParameters, parallelCores, 'none', true, waveformT, targetFilt, [1 1 1 1], [], targetProgress);
% Note: F(2) is only a distinct frequency-spectrum figure when filtering
% is on; otherwise detectMinis returns just the curation figure and F(1)
% and F(2) are the same handle (MATLAB scalar-expands the single output
% into both slots). Previously this was closed here unconditionally,
% which (when filtering was off) silently closed the target curation
% figure itself. It is now left open and handled by the caller, which
% keeps the curation figure and folds any real spectrum figure into the
% combined results panel.
if isempty(targetSpectra)
    targetFilt.state = 'spectrum';
    targetFilt.nSweeps = targetProperties.hd.lActualEpisodes;
    [~, targetSpectra] = filterMinis(targetProperties.sweep, targetProperties.dt, targetFilt, false);
end
if isfield(targetWave,'parameters')
    targetMinis = [targetMinis repmat([targetWave.parameters.averageAmp targetWave.parameters.medianAmp targetWave.parameters.tau_m], size(targetMinis,1),1)];
end
[noiseMinis, ~, ~, noiseProperties, F(3:4), noiseSpectra] = detectMinisStandalone(false, noiseFilename, detectionParameters, noiseExcludedTimes,...
    classificationParameters, parallelCores, 'none', true, waveformN, noiseFilt, [1 1 1 0], [], noiseProgress);
% See the note above F(2): the same applies to F(4) for the noise pass.
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

% Incidence rate (events/s within each 1-second window, retrospectively
% assigned to every event in that window -- see incidenceRate.m):
targetDurationMs = numel(targetProperties.sweep) * targetProperties.dt;
noiseDurationMs = numel(noiseProperties.sweep) * noiseProperties.dt;
targetIncidenceRates = incidenceRate(targetMinis(:,2), targetDurationMs, excludedTimesToIntervals(targetExcludedTimes));
noiseIncidenceRates = incidenceRate(noiseMinis(:,2), noiseDurationMs, excludedTimesToIntervals(noiseExcludedTimes));