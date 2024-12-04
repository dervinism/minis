% Use minisHeadles function to execute 'minis' programatically using
% Matlab. Input and output variables of the function are described below. A
% usage example of how to set up the input variable is given at the end of
% the file.
%
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

input.task = 'detectionHeadless';
input.loadTargetFileInput = 'target_file_path.abf';
input.loadNoiseFileInput = 'noise_file_path.abf';
input.tau_PSPm = '...';
input.tau_m = '13.5288';
input.L = '0.6';
input.loSimAmp = '0.01';
input.maxAmpBottomErr = '726.387';
input.maxAmpBottomDev = '136.34';
input.maxAmpTopErr = '120.003';
input.maxAmpMidErr = '247.006';
input.maxAmpMidDev = '36.3331';
input.maxAmpTopDev = '19.9999';
input.SDupbound = '...';
input.maxErr = '2260.41';
input.SDlobound = '...';
input.maxRTErr = '939.702';
input.maxAmpErr = '1001.71';
input.maxDevAmp = '201.701';
input.maxDevRT = '333.014';
input.maxDev = '196.339';
input.distBaseline = 'Subtracted'; % 'Zero' or 'Subtracted'
input.distType = 'Normal'; % 'Normal', 'Bimodal Normal', 'Trimodal Normal', 'Quadrimodal Normal', 'Log Normal', 'Bimodal Log Normal', 'Trimodal Log Normal', 'Gaussian', 'Bimodal Gaussian', 'Trimodal Gaussian', 'Quadrimodal Gaussian', 'Skew-normal', 'Bimodal Skew-normal', 'Trimodal Skew-normal'
input.voltageClamp = false; % true or false
input.downGoing = false; % true or false
input.pulseDuration = '0.5';
input.RTbinSize = '0.5'; % '0.25' or '0.5'
input.RTint = '10-90%'; % '10-90%' or '20-80%'
input.endGlitchNoise = '...';
input.startGlitchNoise = '...';
input.endPulseNoise = '1.1,2.2';
input.startPulseNoise = '0.05,1.65';
input.endGlitchTarget = '...';
input.startGlitchTarget = '...';
input.endPulseTarget = '1.1,2.2';
input.startPulseTarget = '0.05,1.65';
input.smoothWindow = '1.5';
input.Ampupbound = '10';
input.Amplobound = '0.1';
input.peakIntegrationPeriod = '2.5';
input.baselineDuration = '2';
input.maxTimeToPeak = '10';

input.filtering = 'on';
input.filtfs = '50, 150';

input.options.bounds = [0.3 0.05 500 0.5 2.8 0 0      0      -1 0      0      -0.0001 0      0      -1      0      0      -0.0001 0      0      -1 0      0      -0.0001;
                        0.3 0.05 500 0.5 2.8 0 0.0001 0.0001  1 0.0010 0.0010  0.0001 0.0001 0.0001  1.0000 0.0010 0.0010  0.0001 0.0001 0.0001  1 0.0010 0.0010  0.0001];
input.options.nGenerations = 200;
input.options.parallelCores = 'max';
input.options.fullParallel = 1;
input.options.tauRange = 1;
input.options.cluster = 0;
input.options.cliff = 0;
input.options.figureDisplay = 1;
input.options.clusterProfile = 'local';
input.options.SDlobound = 0;
input.options.SDupbound = 0;

[status, output] = minisHeadless(input);