function input = settings2input(settings)
% input = settings2input(settings)
%
% Function settings2input converts minis settings file to an input
% structure variable used to run minis in a headless mode.
% Input: settings - a path to the settings file.
% Output: input - a structure variable with fields as described in the
%                 documentation of the minisHeadless function.

% Main parameters
if ~isfield(settings, 'task')
  input.task = 'detectionHeadless';
else
  input.task = settings.task;
end
input.loadTargetFileInput = settings.loadTargetFileInput;
input.loadNoiseFileInput = settings.loadNoiseFileInput;
input.tau_PSPm = settings.tau_PSPmEdit;
input.tau_m = settings.tau_mEdit;
input.L = settings.LEdit;
input.loSimAmp = settings.loSimAmpEdit;
input.maxAmpBottomErr = settings.maxAmpBottomErrEdit;
input.maxAmpBottomDev = settings.maxAmpBottomDevEdit;
input.maxAmpTopErr = settings.maxAmpTopErrEdit;
input.maxAmpMidErr = settings.maxAmpMidErrEdit;
input.maxAmpMidDev = settings.maxAmpMidDevEdit;
input.maxAmpTopDev = settings.maxAmpTopDevEdit;
input.SDupbound = settings.SDupbound;
input.maxErr = settings.maxErrEdit;
input.SDlobound = settings.SDlobound;
input.maxRTErr = settings.maxRTErrEdit;
input.maxAmpErr = settings.maxAmpErrEdit;
input.maxDevAmp = settings.maxDevAmpEdit;
input.maxDevRT = settings.maxDevRTEdit;
input.maxDev = settings.maxDevEdit;
if settings.distBaselineEdit == 1
  input.distBaseline = "Zero";
else
  input.distBaseline = "Subtracted";
end
if settings.distTypeEdit == 1
  input.distType = "Normal";
elseif settings.distTypeEdit == 2
  input.distType = "Bimodal Normal";
elseif settings.distTypeEdit == 3
  input.distType = "Trimodal Normal";
elseif settings.distTypeEdit == 4
  input.distType = "Quadrimodal Normal";
elseif settings.distTypeEdit == 5
  input.distType = "Log Normal";
elseif settings.distTypeEdit == 6
  input.distType = "Bimodal Log Normal";
elseif settings.distTypeEdit == 7
  input.distType = "Trimodal Log Normal";
elseif settings.distTypeEdit == 8
  input.distType = "Gaussian";
elseif settings.distTypeEdit == 9
  input.distType = "Bimodal Gaussian";
elseif settings.distTypeEdit == 10
  input.distType = "Trimodal Gaussian";
elseif settings.distTypeEdit == 11
  input.distType = "Quadrimodal Gaussian";
elseif settings.distTypeEdit == 12
  input.distType = "Skew-normal";
elseif settings.distTypeEdit == 13
  input.distType = "Bimodal Skew-normal";
elseif settings.distTypeEdit == 14
  input.distType = "Trimodal Skew-normal";
end
input.voltageClamp = settings.voltageClampCheckbox;
input.downGoing = settings.downGoingCheckbox;
input.pulseDuration = settings.pulseDurationEdit;
if settings.RTbinSizeEdit == 1
  input.RTbinSize = 0.25;
else
  input.RTbinSize = 0.5;
end
if settings.RTintEdit == 1
  input.RTint = '10-90%';
else
  input.RTint = '20-80%';
end
input.endGlitchNoise = settings.endGlitchNoiseEdit;
input.startGlitchNoise = settings.startGlitchNoiseEdit;
input.endPulseNoise = settings.endPulseNoiseEdit;
input.startPulseNoise = settings.startPulseNoiseEdit;
input.endGlitchTarget = settings.endGlitchTargetEdit;
input.startGlitchTarget = settings.startGlitchTargetEdit;
input.endPulseTarget = settings.endPulseTargetEdit;
input.startPulseTarget = settings.startPulseTargetEdit;
input.smoothWindow = settings.smoothWindowEdit;
input.Ampupbound = settings.AmpupboundEdit;
input.Amplobound = settings.AmploboundEdit;
input.peakIntegrationPeriod = settings.peakIntegrationPeriodEdit;
input.baselineDuration = settings.baselineDurationEdit;
input.maxTimeToPeak = settings.maxTimeToPeakEdit;

if ~isfield(settings, 'filtering')
  input.filtering = 'off';
else
  input.filtering = settings.filtering;
end
if ~isfield(settings, 'filtfs')
  input.filtfs = '50, 150';
else
  input.filtfs = settings.filtfs;
end

% Options
input.options = settings.options;

% Extra distrobution fitting parameters
if isfield(settings.options, 'estimateTauLo')
  input.options.estimateTauLo = settings.options.estimateTauLo;
else
  input.options.estimateTauLo = true;
end
if isfield(settings.options, 'estimateTauHi')
  input.options.estimateTauHi = settings.options.estimateTauHi;
else
  input.options.estimateTauHi = true;
end
if isfield(settings.options, 'optimisationData')
  input.options.optimisationData = settings.options.optimisationData;
else
  input.options.optimisationData = '';
end
if isfield(settings.options, 'resumeOptimisation')
  input.options.resumeOptimisation = settings.options.resumeOptimisation;
else
  input.options.resumeOptimisation = '';
end