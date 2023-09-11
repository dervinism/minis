function [status, output] = minisHeadless(input)
% [status, output] = minisHeadless(input)
% minisHeadless function executes minis software in Matlab environment.
%
% Input: input - a structure with the following fields:
%
%         'task' - A task to carry out (character array). Available task
%           options are:
%           'preprocess' - preprocess current-clamp recording data.
%           'detection' - carry out detection of spontaneous postsynaptic
%             potentials/currents in the whole-cell patch-clamp recording
%             data.
%           'detectionHeadless' - carry out detection of spontaneous
%             postsynaptic potentials/currents in the whole-cell
%             patch-clamp recording data but without displaying any
%             graphical windows or saving any of the data. In this case
%             detection results are stored in the 'output' (1) variable.
%             This mode can be useful when embedding minis detection
%             algorithm within custom-written code.
%           'detectCompare' - carry out detection of spontaneous
%             postsynaptic potentials/currents in the whole-cell
%             patch-clamp recording data stored in target and noise files
%             and compare the two.
%           'errorBounds' - estimate error bounds for automated simulated
%             spontaneous postsynaptic potential distribution fitting.
%           'autoDistributionFit' - carry out an automated simulated
%             spontaneous postsynaptic potential distribution fitting.
%           'simulation' - carry out simulation of spontaneous
%             postsynaptic potentials followed by detection and detection
%             performance evaluation.
%           'simulationHeadless' - carry out simulation of spontaneous
%             postsynaptic potentials followed by detection and detection
%             performance evaluation without invoking a GUI and without
%             saving ABFs with simulated traces. In this case the
%             simulation traces and detection performance results are
%             stored in the 'output' (2) variable.
%
%         'loadTargetFileInput' - A full path to a target file (character
%           array). If you are using minis in a detectionHeadless mode, you
%           can supply a series data vector instead. In this case
%           'loadTargetFileInput' should be a structure with the following
%           fields:
%             'dt' - a sampling interval in milliseconds (scalar).
%             'lActualEpisodes' - a number of recordings sweeps (scalar).
%             'sweep' - membrane potential data row vector (mV). You can
%               populate with zeros if running in a voltage clamp mode.
%             'current' - membrane current data row vector (nA). You can
%               populate with zeros if running in a current clamp mode.
%         'loadNoiseFileInput' - A full path to a noise file (character
%           array). If you are using minis in a simulationHeadless mode,
%           you can supply a series data vector instead. In this case
%           'loadNoiseFileInput' should be a structure with the following
%           fields:
%             'dt' - a sampling interval in milliseconds (scalar).
%             'lActualEpisodes' - a number of recordings sweeps (scalar).
%             'sweep' - membrane potential data row vector (mV).
%             'current' - membrane current data row vector (nA). You can
%               also populate with zeros.
%
%         'maxTimeToPeak' - the beginning of the search window for finding
%           the amplitude, and the rise time. It is measured in time to the
%           peak, ms (character array);
%         'baselineDuration' - the default duration of the baseline of a
%           detected event, ms (character array);
%         'peakIntegrationPeriod' - the duration of the
%           refractory/integration period of a minis-like event, ms
%           (character array);
%         'Amplobound' - the lower bound on amplitude of a minis-like
%           event, mV (character array);
%         'Ampupbound' - the upper bound on amplitude of a minis-like
%           event, mV (character array);
%         'smoothWindow' - the element-wise size of the normal smoothing
%           window (character array);
%         'RTint' - a character array representing the rise time interval
%           of choice: '10-90%' for a 10-90% rise time (default) and
%           '20-80%' for a 20-80% rise time.
%         'downGoing' - a logical representing the direction of synaptic
%           events (up- or down-going). If true, events are treated as
%           down-going and the voltage or current trace is inverted prior
%           to running detection.
%         'voltageClamp' - a logical indicating that voltage clamp data is
%           being analysed. If true, detection is performed on the current
%           channel instead of the voltage channel (default is current
%           clamp).
%
%         'RTbinSize' - Rise time classification histogram bin size in ms.
%           It can either be '0.25' or '0.5' (character array; default is
%           0.5).
%
%         'startPulseTarget' - Beginning of (the first) pulse in the target
%           file (s; character array; optional).
%         'endPulseTarget' - End of (the first) pulse in the target file
%           (s; character array; optional).
%         'startGlitchTarget' - End of a glitch period in the target file
%           (s; character array; optional).
%         'endGlitchTarget' - End of a glitch period in the target file (s;
%           character array; optional).
%         'startPulseNoise' - Beginning of (the first) pulse in the noise
%           file (s; character array; optional).
%         'endPulseNoise' - End of (the first) pulse in the noise file (s;
%           character array; optional).
%         'startGlitchNoise' - End of a glitch period in the noise file (s;
%           character array; optional).
%         'endGlitchNoise' - End of a glitch period in the noise file (s;
%           character array; optional).
%         'pulseDuration' - Brief (second) pulse duration (ms; character
%           array). Default value is '0.5'.
%
%         'loSimAmp' - Simulated amplitude lower bound (mV; character
%           array).
%         'L' - the starting size of the electronic length of the dendritic
%           cylinder in simulations (unitless). It corresponds to the real
%           dendritic cylinder length (measured in micrometers) divided by
%           the dendritic length constant (l/λ). The default value is 0.6.
%         'tau_m' - Passive membrane time constant (ms), lower limit
%           (character array).
%         'tau_PSPm' - Passive membrane time constant (ms), upper limit
%           (character array).
%
%         'distType' - Simulated minis distribution type fitted to
%           experimental data during optimisation (character array).
%           Available types are 'Normal', 'Bimodal Normal',
%           'Trimodal Normal', 'Quadrimodal Normal', 'Log Normal',
%           'Bimodal Log Normal', 'Trimodal Log Normal', 'Gaussian',
%           'Bimodal Gaussian', 'Trimodal Gaussian',
%           'Quadrimodal Gaussian', 'Skew-normal', 'Bimodal Skew-normal',
%           'Trimodal Skew-normal'. 'Quadrimodal Normal' is the default
%           option.
%         'distBaseline' - Fitted distribution baseline (character array).
%           Options are either 'Zero' (default) or 'Subtracted'. This
%           option allows the user to add to the simulated distribution all
%           the positive events after subtracting the two-dimensional
%           amplitude and rise time noise distribution from the target
%           distribution (Subtracted option). The default is Zero which
%           does not add anything.
%         'SDupbound' - Membrane potential standard deviation upper bound
%           based on 15 ms size window average (mV; character array). This
%           is optional and if set to be an empty character array, it is
%           not used in optimisation (by default).
%         'SDlobound' - Membrane potential standard deviation lower bound
%           based on 15 ms size window average (mV; character array). This
%           is optional and if set to be an empty character array, it is
%           not used in optimisation (by default).
%         'maxErr' - Maximum combined SAD (character array).
%         'maxAmpErr' - Maximum amplitude SAD (character array).
%         'maxRTErr' - Maximum rise time SAD (character array).
%         'maxAmpBottomErr' - Maximum amplitude top 50% SAD (character
%           array).
%         'maxAmpMidErr' - Maximum amplitude top 10% SAD (character array).
%         'maxAmpTopErr' - Maximum amplitude top 2% SAD (character array).
%         'maxDev' - Maximum combined MAD (character array).
%         'maxDevAmp' - Maximum amplitude MAD (character array).
%         'maxDevRT' - Maximum rise time MAD (character array).
%         'maxAmpBottomDev' - Maximum amplitude top 50% MAD (character
%           array).
%         'maxAmpMidDev' - Maximum amplitude top 10% MAD (character array).
%         'maxAmpTopDev' - Maximum amplitude top 2% MAD (character array).
%
%         'filtering' - a character array that can be either set to 'on' or
%           'off' (default). If set to 'on', the electrophysiological data
%           is band-stop filtered to remove frequency components defined in
%           'filtfs'. This parameter is used only when running minis in a
%           detectionHeadless or simulationHeadless modes.
%         'filtfs' - frequencies to be filtered out if 'filtering' is set
%           to 'on' (optional). It is a cell with a character array inside
%           the brackets. The default value is {'50, 150'}. This parameter
%           is used only when running minis in a detectionHeadless or
%           simulationHeadless modes.
%
%         'options' - A structure variable to store optimisation options.
%           It contains the following fields:
%           'bounds' - a row matrix with the lower (top row) and upper
%             (bottom row) bounds on the simulated distribution parameters
%             used during automated simualted distribution fitting.
%           'nGenerations' - a number of genetic algorithm generations
%             (scalar).
%           'parallelCores' - task parallelisation (number of cores;
%             character array). When entering the number of parallel cores,
%             instead of giving a number, the user can enter 'min' or 'max'
%             which would give the minimum or maximum available cores on
%             the cluster profile (computer), respectively.
%           'fullParallel' - a logical with true corresponding to the
%             maximum available number of parallel cores and false
%             (default) using the number of cores specified in the
%             'parallelCores' field.
%           'tauRange' - a logical with true corresponding having the
%             passive membrane time constant as an extra parameter
%             controlled by the optimisation procedure. Default is false.
%           'cluster' - true or false for evaluating on a remote parallel
%             cluster (false by default).
%           'clusterProfile' - specify the parallel cluster profile name
%             (character array). Default is 'local'.
%           'cliff' - a logical with true corresponding to using cliff
%             constraint (refer to 'minis' documentation; false by default).
%           'figureDisplay' - a logical for displaying figures during the
%             optimisation process (true by default).
%           'SDlobound' - a logical for using lower bound on membrane
%             potential standard deviation as another controlled parameter
%             during optimisation.
%           'SDupbound' - a logical for using upper bound on membrane
%             potential standard deviation as another controlled parameter
%             during optimisation.
%           'estimateTauLo' - a logical for estimating passive membrane
%             time constant lower limit based on measured decays of
%             detected minis-like events. This field is used only during
%             autoDistributionFitHeadless task. If specified to be true,
%             the supplied range of membrane time constants in the
%             Simulation Paprameters panel is overridden and a new range is
%             established if 'tauRange' was set to be true.
%           'estimateTauHi' - a logical for estimating passive membrane
%             time constant upper limit based on measured decays of
%             detected minis-like events. This field is used only during
%             autoDistributionFitHeadless task. If specified to be true,
%             the supplied range of membrane time constants in the
%             Simulation Paprameters panel is overridden and a new range is
%             established if 'tauRange' was set to be true.
%           'optimisationData' - a character array containing optimisation
%             data file path obtained in error bound estimation. Supply
%             only when running autoDistributionFitHeadless task. It is an
%             empty character array by default.
%           'resumeOptimisation' - a character array containing unfinished
%             optimisation data file path obtained during previous
%             optimisation. Supply only when running
%             autoDistributionFitHeadless task. It is an empty character
%             array by default.
% Output: output (1) - a structure with the following fields (produced only
%           during the detectionHeadless mode, otherwise empty):
%           'minisArray' - a matrix with rows corresponding to the
%             detected minis-like events sorted according their temporal
%             order and columns (from left to right) corresponding to the
%             following characteristics of detected event:
%               (1)  - the peak value, mV or nA;
%               (2)  - the peak time, ms;
%               (3)  - the peak index;
%               (4)  - the amplitude, mV;
%               (5)  - the baseline value, mV;
%               (6)  - the start of a baseline, ms;
%               (7)  - the end of a baseline, ms;
%               (8)  - the element-wise start of a baseline, ms;
%               (9)  - the element-wise end of a baseline, ms;
%               (10) - the element-wise 0-100% rise time (from the start of
%                      the baseline to the peak).
%               (11) - the 0-100% rise time, ms.
%               (12) - the 10-90% rise time, ms.
%               (13) - the 10% rise time mark, ms.
%               (14) - the 50% rise time mark, ms.
%               (15) - the 90% rise time mark, ms.
%               (16) - the index of the 10% rise time mark.
%               (17) - the index of the 50% rise time mark.
%               (18) - the index of the 90% rise time mark.
%               (19) - membrane potential or current value at the 10% rise
%                      time mark, mV or nA.
%               (20) - membrane potential or current value at the 50% rise
%                      time mark, mV or nA.
%               (21) - membrane potential or current value at the 90% rise
%                      time mark, mV or nA.
%               (22) - 1/e decay, ms.
%             'waveform' - a structure containing the fields:
%               'averageTrace' is an average recording trace in mV or nA.
%               'riseTimeArray' is a vector with rise times corresponding
%                 to counts stored in the 'riseTimeDist' vector.
%               'riseTimeDist' is a vector with a rise time distribution of
%                 large (top 10%) amplitude minis-like events that were
%                 used for obtaining the average waveform.
%               'parameters' is a structure variable containing estimated
%                 waveform parameters 'peak' (mV or nA), 'BL' (baseline, mV
%                 or nA), 'Amp' (amplitude, mV or nA), 'risetime' (ms),
%                 'tau_m' (effective dendritic membrane time constant, ms).
%                 Additional fields like averageAmp and medianAmp are also
%                 added which correspond to the mean and median amplitude
%                 of the top 10% of the largest detected minis-like events.
%             'spectrum' - a matrix composed of four row vectors. The first
%               vector contains amplitudes. The second one contains a
%               power spectrum which should correspond to decibels (dB).
%               The third vector contains phase information. The fourth
%               vector contains frequencies.
%         output (2) - a structure with the following fields (produced only
%           during the simulationHeadless mode, otherwise empty):
%           'simData' - a row matrix with the top row being the noise +
%             simulated membrane potential data after band-stop filtering
%             and smoothing performed by the detection algorithms and the
%             bottom row being the current data.
%           'simDataRaw' - a row matrix with the top row being the noise +
%             simulated membrane potential data after band-stop filtering
%             and the bottom row being the current data.
%           'detectionPerformance' - a structure variable with the
%             following fields:
%             'sensitivity' - true positive (hit) rate.
%             'specificity' - specificity or 1 - false positive (false
%               alarm) rate.
%             'FPR' - false positive rate.
%             'dPrime' - sensitivity index.
%             'performance' - a row matrix containing logical indices
%               corresponding to
%                 row 1: simulated event positions (peaks);
%                 row 2: hits (detected positions) + misses;
%                 row 3: hits (detected positions) + false alarms;
%                 row 4: hits (detected positions);
%                 row 5: misses;
%                 row 6: false alarms;
%                 row 7: correct rejections.
%             'falseI' - locations (indices) of prominent noise events.
%             'falseT' - times of prominent noise events
%           'detectionParameters' - a structure variable containing part
%             of input parameters controling the detection task.
%           'simulationParameters' - a structure variable containing part
%             of input parameters controling simulation of sEPSPs.
%           'optimisationParameters' - a structure variable containing
%             part of input parameters controling the automated
%             distribution fitting task.
%           'classificationParameters' - a structure variable containing
%             part of input parameters controling distribution binning.
%           'filtfs' - band-stop filtering stop frequencies.
%           'simulatedEventInfo' - a matrix with rows corresponding to
%             simulated events and columns corresponding to (1) minis
%             count, (2) amplitude, (3) rise time, (4) electrotonic charge
%             input site distance from the measuring site, (5) electrotonic
%             length of the dendritic cylinder, (6) onset index, (7) onset
%             time, (8) peak index, (9) peak time.
arguments
  input struct
end
                      
% Initialise output variables
status = -1;
output = [];

% Validate input arguments
input = validateInput(input);

% Save the working directory
if (~isdeployed)
  currentFullFilename = mfilename('fullpath');
  currentPath = fileparts(currentFullFilename);
  addpath(genpath(currentPath));
else
  if ispc
    [~, currentPath] = system('cd');
  elseif isunix
    [~, currentPath] = system('pwd');
  end
end
input.wd = currentPath;
input.ld = currentPath;

% Initialise general variables:
graphicsFormats = {'*.fig','Matlab figure (*.fig)'
    '*.ai','Adobe Illustrator `88 (*.ai)'
    '*.bmp','Windows bitmap (*.bmp)'
    '*.emf','Enhanced metafile (*.emf)'
    '*.eps','EPS Level 1 (*.eps)'
    '*.jpg','JPEG image (*.jpg)'
    '*.m','MATLAB file (*.m)'
    '*.pbm','Portable bitmap (*.pbm)'
    '*.pcx','Paintbrush 24-bit (*.pcx)'
    '*.pdf','Portable Document Format (*.pdf)'
    '*.pgm','Portable Graymap (*.pgm)'
    '*.png','Portable Network Graphics (*.png)'
    '*.ppm','Portable Pixmap (*.ppm)'
    '*.tif','TIFF image, compressed (*.tif)'};

[initialised, detectionParameters] = initDetectParam(input);
if ~initialised
    return
end

classificationParameters = initClassParam(input, detectionParameters.Ampupbound);

parallelCores = initParallelisation(input.options.parallelCores);

% Determine the choice of action:
task = input.task;


if strcmpi(task, 'preprocess')
    % Pulses:
    [initialised, startPulse, endPulse] = initPulseTarget(input);
    if ~initialised
        return
    end
    
    % Pre-process data:
    preprocessMinis(startPulse, endPulse, detectionParameters, classificationParameters, input.ld, graphicsFormats, parallelCores);
    
elseif strcmpi(task, 'detection')
    % Target file:
    targetFilename = input.loadTargetFileInput;
    if strcmpi(targetFilename, '>>> <<<') || isempty(targetFilename)
        errmsgNoFile = 'Error: no target file supplied';
        msgbox(errmsgNoFile,'Error','Error');
        return
    end
    
    % Pulses and glitches:
    [initialised, targetExcludedTimes] = initExclTimesTarget(input);
    if ~initialised
        return
    end
    
    % Estimate tau_m based on impulses for later using it as a reference time interval for tau_m estimation based on spontaneous PSPs:
    if detectionParameters.voltageClamp
        waveform.estimate = 1;
        waveform.riseTimeArrayExt = classificationParameters.riseTimeArray;
        waveform.tau_m = 10;
    else
        waveform = initWaveform(input, detectionParameters.pulseDuration, classificationParameters.riseTimeArray, targetExcludedTimes, targetFilename);
    end
    
    % Filtering:
    filtering = noiseFilterDlg;
    filtering.excludedTimes = targetExcludedTimes;
    
    % Detect minis-like events in a file of choice:
    [minisArray, ~, waveform, ~, F] = detectMinisStandalone(true, targetFilename, detectionParameters, targetExcludedTimes, classificationParameters,...
        parallelCores, 'none', true, waveform, filtering, [1 1 1 1]);
    if ~isempty(waveform)
        minisArray = [minisArray repmat([waveform.parameters.averageAmp waveform.parameters.medianAmp waveform.parameters.tau_m], size(minisArray,1), 1)];
    end
    
    % Save the detected event log and figures:
    input.ld = saveTargetEventLog(minisArray, F, waveform, filtering, detectionParameters.RTinterval, graphicsFormats, input.ld);
    
elseif strcmpi(task, 'detectionHeadless')
    % Target file:
    targetFilename = input.loadTargetFileInput;
    if strcmpi(targetFilename, '>>> <<<') || isempty(targetFilename)
        errmsgNoFile = 'Error: no target file supplied';
        msgbox(errmsgNoFile,'Error','Error');
        return
    end
    
    % Pulses and glitches:
    [initialised, targetExcludedTimes] = initExclTimesTarget(input);
    if ~initialised
        return
    end
    
    % Estimate tau_m based on impulses for later using it as a reference time interval for tau_m estimation based on spontaneous PSPs:
    if detectionParameters.voltageClamp
        waveform.estimate = 2;
        waveform.riseTimeArrayExt = classificationParameters.riseTimeArray;
        waveform.tau_m = 10;
    else
        waveform = initWaveform(input, detectionParameters.pulseDuration, classificationParameters.riseTimeArray, targetExcludedTimes, targetFilename);
        waveform.estimate = 2;
    end
    
    % Filtering:
    filtering.state = input.filtering;
    filtering.excludedTimes = targetExcludedTimes;
    filtering.filtfs = input.filtfs;
    
    % Detect minis-like events in a file of choice:
    [minisArray, ~, waveform, ~, ~, spectrum] = detectMinisStandalone(false, targetFilename, detectionParameters, targetExcludedTimes, classificationParameters,...
        parallelCores, 'none', false, waveform, filtering, [1 1 1 1], false);
    if ~isempty(waveform)
        minisArray = [minisArray repmat([waveform.parameters.averageAmp waveform.parameters.medianAmp waveform.parameters.tau_m], size(minisArray,1), 1)];
    end
    
    % Assign output:
    output.minisArray = minisArray;
    output.waveform = waveform;
    output.waveform.riseTimeArray = output.waveform.riseTimeArrayExt(2:end);
    output.waveform = rmfield(output.waveform, 'riseTimeArrayExt');
    output.spectrum = spectrum;
    
elseif strcmpi(task, 'detectCompare')
    % Target file:
    targetFilename = input.loadTargetFileInput;
    if strcmpi(targetFilename, '>>> <<<') || isempty(targetFilename)
        errmsgNoFile = 'Error: no target file supplied';
        msgbox(errmsgNoFile,'Error','Error');
        return
    end
    
    % Noise file:
    noiseFilename = input.loadNoiseFileInput;
    if strcmpi(noiseFilename, '>>> <<<') || isempty(noiseFilename)
        errmsgNoFile = 'Error: no noise file supplied';
        msgbox(errmsgNoFile,'Error','Error');
        return
    end
    
    % Pulses and glitches:
    [initialised, targetExcludedTimes] = initExclTimesTarget(input);
    if ~initialised
        return
    end
    
    [initialised, noiseExcludedTimes] = initExclTimesNoise(input);
    if ~initialised
        return
    end
    
    % Estimate tau_m based on impulses for later using it as a reference time interval for tau_m estimation based on spontaneous PSPs:
    waveform = initWaveform(input, detectionParameters.pulseDuration, classificationParameters.riseTimeArray, targetExcludedTimes, targetFilename);
    
    % Filtering:
    filtering = noiseFilterDlg;
    
    % Detect minis-like events in both target and noise files and compare them:
    [minisTarget, minisNoise, waveform, F] = compareMinis(targetFilename, noiseFilename, targetExcludedTimes, noiseExcludedTimes, detectionParameters,...
        classificationParameters, filtering, waveform, parallelCores);
    
    % Save the detected event logs and figures:
    input.ld = saveEventLog(minisTarget, minisNoise, F, waveform, filtering, detectionParameters.RTinterval, graphicsFormats, input.ld);
    
elseif strcmpi(task, 'errorBounds')
    % Pulses:
    [initialised, startPulse, endPulse] = initPulseTarget(input);
    if ~initialised
        return
    end
    targetExcludedTimes = struct('startPulse', startPulse, 'endPulse', endPulse);
    
    % Estimate error margin for optimisation:
    errorMinis(targetExcludedTimes, detectionParameters, classificationParameters, input.ld, input.wd, parallelCores);
    
elseif strcmpi(task, 'autoDistributionFit')
    % Target file:
    targetFilename = input.loadTargetFileInput;
    if strcmpi(targetFilename, '>>> <<<') || isempty(targetFilename)
        errmsgNoFile = 'Error: no target file supplied';
        msgbox(errmsgNoFile,'Error','Error');
        return
    end
    
    % Noise file:
    noiseFilename = input.loadNoiseFileInput;
    if strcmpi(noiseFilename, '>>> <<<') || isempty(noiseFilename)
        errmsgNoFile = 'Error: no noise file supplied';
        msgbox(errmsgNoFile,'Error','Error');
        return
    end
    
    % Pulses and glitches:
    [initialised, targetExcludedTimes] = initExclTimesTarget(input);
    if ~initialised
        return
    end
    
    [initialised, noiseExcludedTimes] = initExclTimesNoise(input);
    if ~initialised
        return
    end
    
    % Optimisation parameters:
    [initialised, optimisationParameters] = initOptimParam(input);
    if ~initialised
        return
    end
    
    % Simulation parameters:
    [initialised, simulationParameters] = initSimParam(input, detectionParameters.pulseDuration, classificationParameters.riseTimeArray,...
        targetExcludedTimes, targetFilename);
    if ~initialised
        return
    end
    
    % Filtering:
    filtering = noiseFilterDlg;
    filtering.targetExcludedTimes = targetExcludedTimes;
    filtering.noiseExcludedTimes = noiseExcludedTimes;
    
    % Start automatic distribution fitting:
    optimiseMinis(input.wd, targetFilename, noiseFilename, targetExcludedTimes, noiseExcludedTimes, detectionParameters, simulationParameters,...
        optimisationParameters, classificationParameters, filtering, 'TwoDs', parallelCores);
elseif strcmpi(task, 'autoDistributionFitHeadless')
    % Target file:
    targetFilename = input.loadTargetFileInput;
    if strcmpi(targetFilename, '>>> <<<') || isempty(targetFilename)
        errmsgNoFile = 'Error: no target file supplied';
        msgbox(errmsgNoFile,'Error','Error');
        return
    end
    
    % Noise file:
    noiseFilename = input.loadNoiseFileInput;
    if strcmpi(noiseFilename, '>>> <<<') || isempty(noiseFilename)
        errmsgNoFile = 'Error: no noise file supplied';
        msgbox(errmsgNoFile,'Error','Error');
        return
    end
    
    % Pulses and glitches:
    [initialised, targetExcludedTimes] = initExclTimesTarget(input);
    if ~initialised
        return
    end
    
    [initialised, noiseExcludedTimes] = initExclTimesNoise(input);
    if ~initialised
        return
    end
    
    % Optimisation parameters:
    [initialised, optimisationParameters] = initOptimParam(input);
    if ~initialised
        return
    end
    
    % Simulation parameters:
    [initialised, simulationParameters] = initSimParam(input, detectionParameters.pulseDuration, classificationParameters.riseTimeArray,...
        targetExcludedTimes, targetFilename);
    if ~initialised
        return
    end
    
    % Filtering:
    filtering.state = input.filtering;
    filtering.targetExcludedTimes = targetExcludedTimes;
    filtering.noiseExcludedTimes = noiseExcludedTimes;
    filtering.filtfs = input.filtfs;
    
    % Start automatic distribution fitting:
    optimisationParameters.options.figureDisplay = false;
    optimisationParameters.options.headless = true;
    optimiseMinis(input.wd, targetFilename, noiseFilename, targetExcludedTimes, noiseExcludedTimes, detectionParameters, simulationParameters,...
        optimisationParameters, classificationParameters, filtering, 'TwoDs', parallelCores);
elseif strcmpi(task, 'simulation')
    % Noise file:
    noiseFilename = input.loadNoiseFileInput;
    if strcmpi(noiseFilename, '>>> <<<') || strcmpi(noiseFilename, '...') || strcmpi(noiseFilename, '') || isempty(noiseFilename)
        noiseFilename = [];
    end
    
    % Pulses and glitches:
    [initialised, noiseExcludedTimes] = initExclTimesNoise(input);
    if ~initialised
        return
    end
    
    % Optimisation parameters:
    [initialised, optimisationParameters] = initOptimParam(input);
    if ~initialised
        return
    end
    
    % Simulation parameters:
    [initialised, simulationParameters] = initSimParamReduced(input);
    if ~initialised
        return
    end
    
    % Filtering:
    if ~isempty(noiseFilename)
        filtering = noiseFilterDlg;
        filtering.noiseExcludedTimes = noiseExcludedTimes;
    else
        filtering.state = 'off';
        filtering.noiseExcludedTimes = [];
    end
    
    % Simulate minis:
    simulateDetectEvaluate(noiseFilename, noiseExcludedTimes, detectionParameters, simulationParameters, optimisationParameters,...
        classificationParameters, filtering, parallelCores);
      
 elseif strcmpi(task, 'simulationHeadless')
    % Noise file:
    noiseFilename = input.loadNoiseFileInput;
    if strcmpi(noiseFilename, '>>> <<<') || strcmpi(noiseFilename, '...') || strcmpi(noiseFilename, '') || isempty(noiseFilename)
        noiseFilename = [];
    end
    
    % Pulses and glitches:
    [initialised, noiseExcludedTimes] = initExclTimesNoise(input);
    if ~initialised
        return
    end
    
    % Optimisation parameters:
    [initialised, optimisationParameters] = initOptimParam(input);
    if ~initialised
        return
    end
    
    % Simulation parameters:
    [initialised, simulationParameters] = initSimParamReduced(input);
    if ~initialised
        return
    end
    
    % Filtering:
    if ~isempty(noiseFilename)
        filtering.state = input.filtering;
        filtering.excludedTimes = noiseExcludedTimes;
        filtering.filtfs = input.filtfs;
    else
        filtering.state = 'off';
        filtering.noiseExcludedTimes = [];
    end
    
    % Simulate minis:
    output = simulateDetectEvaluate(noiseFilename, noiseExcludedTimes, detectionParameters, simulationParameters, optimisationParameters,...
        classificationParameters, filtering, parallelCores, []);
end
status = 0;