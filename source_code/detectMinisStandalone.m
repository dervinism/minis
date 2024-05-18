function [minis, filterV, waveform, dataProperties, F, spectrum] = detectMinisStandalone(display, varargin)
% DETECTMINISSTANDALONE detects miniature excitatory postsynaptic
% potentials (mEPSPs) or minis (or minis-like events) in the user supplied
% electrophysiological recording data. The recording data could be membrane
% potential in milivolts (mV) or current clamp in nanoamperes (nA).
%
% This function can be used on its own or in conjunction with the minis
% graphical user interface. If used as a standalone function, user should
% type in the display mode, detection parameters, number of recording
% sweeps, excluded times, and the number of parallel processor cores
% intended to use.
%
%   MINIS = DETECTMINISSTANDALONE(display, filename,...
%       detectionParameters, excludedTimes, classificationParameters,...
%       parallelCores, exclTimesAction, edit, waveform, filtering, SD,
%       display2)
%   detects the minis-like events and outputs detection data in MINIS. This
%   variable can later be saved as a text log file (e.g.).DISPLAY is a
%   logical variable with values true for displaying event histograms and
%   false otherwise (default value). FILENAME is the full file name
%   (including a path) of an abf file containing data. DETECTIONPARAMETERS
%   is a structure variable used by the minis-like event detection
%   algorithm with the following fields:
%       'SWstart' - the beginning of the search window for finding the
%           amplitude, and the rise time. It is measured in time to the
%           peak, ms;
%       'BLduration' - the default duration of the baseline of a detected
%           event, ms;
%       'refractoryPeriod' - the duration of the refractory/integration
%           period of a minis-like event, ms;
%       'Amplobound' - the lower bound on amplitude of a minis-like event,
%           mV;
%       'Ampupbound' - the upper bound on amplitude of a minis-like event,
%           mV;
%       'smoothWindow' - the element-wise size of the normal smoothing
%           window;
%       'smoothWindowLite' - the element-wise size of the light smoothing
%           window.
%       'RTinterval' - a string representing the rise time interval of
%           choice: '10-90%' for a 10-90% rise time and '20-80%' for a
%           20-80% rise time.
%       'downGoing' - a logical representing the direction of synaptic
%                     events (up- or down-going). If true, events are
%                     treated as down-going and the voltage or current
%                     trace is inverted prior to running detection.
%       'voltageClamp' - a logical indicating voltage clamp. If true,
%                        detection is performed on the current channel
%                        instead of the voltage channel (default current
%                        clamp).
%   EXCLUDEDTIMES is a vector with times that should not be analysed.
%   CLASSIFICATIONPARAMETERS is a structure variable containing
%   classification ranges for amplitudes and rise times. The fields of this
%   variable are
%       'amplitudeArray' - a vector of classification amplitudes;
%       'amplitudeArrayExt' - a vector of classification amplitudes with an
%           additional bin at the end;
%       'riseTimeArray' - a vector of classification rise times;
%       'riseTimeArrayExt' - a vector of classification rise times with an
%           additional bin at the end;
%   PARALLELCORES is a scalar variable with values set to correspond to the
%   number of processor cores intended to use for the parallel processing
%   in the detection task. EXCLTIMESACTION is a an optional string variable
%   with a value 'crop' for cropping out excluded times. EDIT is an
%   optional variable for controlling whether detected events should be
%   manually edited. A logical true value means that events should be
%   edited and a logical false value is otherwise. The default mode is true.
%   WAVEFORM is a structure variable with fields:
%       'estimate' - a scalar variable that can take values 1, 2, 3, or 0.
%           1 is for estimating the average waveform with the user being
%           prompted for a choice to perform the average waveform
%           estimation, 2 is for the average waveform estimation without
%           prompting the user, 3 is the same as 2 but with a graphical
%           output of estimation, whereas 0 is for not performing the
%           estimation. The default value is 0.
%       'riseTimeArrayExt' - a vector with the rise time range for
%           classification (extended by one element).
%       'tau_m' - a scalar variable containing an estimate of the dendritic
%           membrane time constant based on impulses.
%   FILTERING is a structure variable that is used by the filterMinis
%   function. It contains three fields:
%       'state' - a string that controls the band-stop filtering of the
%           electrophysiological data with values 'on' and 'off'
%           corresponding to filtering and leaving the data unfiltered,
%           respectively. The third option 'spectrum' returns frequency
%           spectrum vector without filtering the data.
%       'pulseEnd' - a scalar corresponding to the end time of the last
%           pulse.
%       'filtfs' - frequencies to be filtered out (optional). It is a cell
%           array of the following example form: {'50, 150'}.
%   SD is a vector variable controlling the estimation of the standard
%   deviation (SD) of the recording trace. The vector contains four
%   elements each of which could be set to either 1 (true) or 0 (false).
%   The following combinations are used to calculate different types of SD:
%       [1 0 0 0] - calculates the SD of the entire recording trace without
%           any adjustments for non-stationarity of the data.
%       [0 1 0 0] - calculates the SD of the entire recording trace with
%           non-stationarity minimised. The calculation is based on
%           15ms-window average.
%       [0 0 1 0] - calculates the SD of the entire recording trace with
%           non-stationarity minimised. The calculation is based on
%           30ms-window average.
%       [0 0 0 1] - calculates the SD of the intervals containing no
%           detected minis-like events. This calculation estimates the SD
%           of the pseudo noise and is based on 15ms-window average.
%       Other combinations are possible allowing to calculate more than one
%           type of SD at the same time.
%   The default combination is [0 0 0 0] where no SD calculations are
%   carried out. If the user specifies any of the combinations containing
%   ones, the estimates are appended to the MINIS output array as
%   additional collumns 22-25 in the order they are listed.
%   DISPLAY2 is a logical controlling the display of figures. If set to
%   false, no figures would be displayed (this does not affect amplitude
%   and rise-time histogram figures). The default is true.
%
%   [MINIS FILTERV] = DETECTMINISSTANDALONE(...)
%   In addition returns the electrophysiological recording vector FILTERV.
%   This option is extremely useful when minis detection is performed as
%   part of the manual or automatic distribution fitting in the minis
%   program. The output vector FILTERV is already filtered and can later be
%   re-used without filtering.
%
%   [MINIS FILTERV WAVEFORM] = DETECTMINISSTANDALONE(...)
%   in addition outputs a structure variables WAVEFORM containing fields:
%       'averageTrace' is an average recording trace in mV or nA.
%       'riseTimeDist' is a vector with a rise time distribution of large
%           amplitude minis-like events that were used for obtaining the
%           average waveform.
%       'parameters' is a structure variable containing estimated waveform
%           parameters 'peak' (mV or nA), 'BL' (baseline, mV or nA), 'Amp'
%           (amplitude, mV or nA), 'risetime' (ms), 'tau_m' (effective
%           dendritic membrane time constant, ms). In case when input
%           variable field waveform.estimate is set to 3, additional fields
%           averageAmp and medianAmp are added which correspond to the mean
%           and median amplitude of the top 10% of the largest detected
%           minis-like events.
%       'F' - handles of the average trace and the rise time distribution
%           figures.
%   If the input waveform variable is not set to 1, the output waveform
%   variable is empty.
%
%   [MINIS FILTERV WAVEFORM DATAPROPERTIES] = DETECTMINISSTANDALONE(...)
%   in addition outputs the recording data obtained from the abf data file.
%   It is intended to use in simulation-optimisation tasks where the noise
%   data needs to be cropped to remove glitches and pulses from the
%   recording trace. In this case DATAPROPERTIES is output already
%   preprocessed.
%
%   [MINIS FILTERV WAVEFORM DATAPROPERTIES F] = DETECTMINISSTANDALONE(...)
%   also outputs a handle vector F of the produced figures.
%
%   [MINIS FILTERV WAVEFORM DATAPROPERTIES F SPECTRUM]
%       = DETECTMINISSTANDALONE(...)
%   In addition outputs the data amplitude, power, and phase spectra.
%   SPECTRUM is a matrix composed of four row vectors. The first vector
%   contains amplitudes. The second vector contains a power spectrum which
%   should correspond to decibels (dB). The third vector contains phase
%   information. The fourth vector contains frequencies.
%



if nargin
    filename = varargin{1};
    if ischar(filename)
        filenameSplit = regexp(filename, '\\*', 'split');
        filenameShort = char(filenameSplit(end));
    else
        filenameShort = filename;
    end
    
    % Load an abf data file:
    if ischar(filename)
        dataProperties = loadABF(filename);
    else
        dataProperties = filename;
        dataProperties.hd.lActualEpisodes = dataProperties.lActualEpisodes;
    end
    
    detectionParameters = varargin{2};
    detectionParameters.sampleInterval = dataProperties.dt;
    detectionParameters.smoothWindow = round(detectionParameters.smoothWindow/detectionParameters.sampleInterval);
    detectionParameters.smoothWindowLite = 8;
    if nargin == 12
        options.display = varargin{11};
    else
        options.display = true;
    end
    if options.display
      options.summaryPlot = true;
    else
      options.summaryPlot = false;
    end
    options.filename = filenameShort;
    
    
    % Determine excluded times:
    excludedTimes = varargin{3};
    if nargin >= 7
        exclTimesAction = varargin{6};
    end
    try
      numSweeps = dataProperties.hd.lActualEpisodes;
    catch
      numSweeps = dataProperties.lActualEpisodes;
    end
    sweepDuration = (length(dataProperties.sweep)*dataProperties.dt - dataProperties.dt)/numSweeps;
    startPulse = excludedTimes.startPulse;
    endPulse = excludedTimes.endPulse;
    startGlitch = excludedTimes.startGlitch;
    endGlitch = excludedTimes.endGlitch;
    if exist('exclTimesAction', 'var') && strcmpi(exclTimesAction, 'crop')
        [dataProperties.sweep, dataProperties.excludedTimes] = cropData(dataProperties.sweep, sweepDuration, numSweeps, 1000*startPulse, 1000*endPulse,...
            1000*startGlitch, 1000*endGlitch, dataProperties.dt);
        [dataProperties.current, dataProperties.excludedTimes] = cropData(dataProperties.current, sweepDuration, numSweeps, 1000*startPulse, 1000*endPulse,...
            1000*startGlitch, 1000*endGlitch, dataProperties.dt);
        excludedTimes = dataProperties.excludedTimes;
    else
        dataProperties.excludedTimes = calcExcludedTimes(sweepDuration, numSweeps, 1000*startPulse, 1000*endPulse, 1000*startGlitch, 1000*endGlitch, dataProperties.dt);
        excludedTimes = dataProperties.excludedTimes;
    end
    
    
    % Determine classification ranges:
    classificationParameters = varargin{4};
    amplitudeArray = classificationParameters.amplitudeArray;
    ampLim2D = (amplitudeArray(2) - amplitudeArray(1))/2;
    riseTimeArray = classificationParameters.riseTimeArray;
    RTLim2D = (riseTimeArray(2) - riseTimeArray(1))/2;
    
    
    % Determine the rest of the input variables:
    parallelCores = varargin{5};
    if nargin >= 8
        options.edit = varargin{7};
        if nargin >= 9
            waveform = varargin{8};
            waveform.nSweeps = numSweeps;
            if nargin >= 10
                filtering = varargin{9};
                filtering.nSweeps = numSweeps;
                if nargin >= 11
                    options.SD = varargin{10};
                end
            else
                filtering.state = 'on';
            end
        else
            waveform.estimate = false;
        end
    else
        options.edit = true;
    end
    clear endGlitch endPulse exclTimesAction filename filenameShort filenameSplit startGlitch startPulse sweepDuration varargin
else
    %% Load an abf data file:
    dataProperties = loadABF;
    
    
    %% Display the histograms of detected events?
    display = false;
    
    
    %% Detection Parameters:
    % Each of these variables is described inside the detectMinis function.
    % Type in the command line "help detectMinis" for info.
    numSweeps = 1;
    options.interactMode = false;
    options.pauseTime = 0.5;
    options.summaryPlot = true;
    options.edit = true;
    options.display = true;
    if ~exist('detectionParameters','var')
        detectionParameters = struct('sampleInterval', dataProperties.dt, 'SWstart', 14, 'BLduration', 2, 'refractoryPeriod', 3,...
            'Amplobound', .02, 'Ampupbound', 1, 'smoothWindow', 30, 'smoothWindowLite', 8, 'nSweeps', numSweeps, 'RTinterval', '10-90%', 'downGoing', false, 'voltageClamp', false);
    end
    
    
    %% Determine excluded times:
    % For the explanation of these variables please see the comments inside the
    % subfunction calcExcludedTimes.
    sweepDuration = (length(dataProperties.sweep)*dataProperties.dt - dataProperties.dt)/numSweeps;
    startPulse = [];
    endPulse = [];
    startGlitch = [];
    endGlitch = [];
    exclTimesAction = 'none';
    if exist('exclTimesAction', 'var') && strcmpi(exclTimesAction, 'crop')
        [dataProperties.sweep, dataProperties.excludedTimes] = cropData(sweepDuration, numSweeps, 1000*startPulse, 1000*endPulse, 1000*startGlitch, 1000*endGlitch, dataProperties.dt);
    else
        excludedTimes = calcExcludedTimes(sweepDuration, numSweeps, 1000*startPulse, 1000*endPulse, 1000*startGlitch, 1000*endGlitch, dataProperties.dt);
    end
    
    
    %% Determine classification ranges:
    ampStepSize = .01;
    RTstepSize = .5;
    amplitudeArray = Amplobound: ampStepSize :.34;
    amplitudeArrayExt = [amplitudeArray amplitudeArray(end)+ampStepSize];
    ampLim2D = (amplitudeArray(2) - amplitudeArray(1))/2;
    riseTimeArray = RTstepSize: RTstepSize :10;
    riseTimeArrayExt = [riseTimeArray riseTimeArray(end)+RTstepSize];
    RTLim2D = (riseTimeArray(2) - riseTimeArray(1))/2;
    
    
    %% Determine the number of parallel processors:
    parallelCores = 1;
    if parallelCores > 1
        parallelCores = floor(parallelCores);
        delete(gcp('nocreate'));
        parpool(parallelCores);
    end
    
    classificationParameters = struct('amplitudeArray', amplitudeArray, 'amplitudeArrayExt', amplitudeArrayExt,...
        'riseTimeArray', riseTimeArray, 'riseTimeArrayExt', riseTimeArrayExt);
    waveform.estimate = 1;
    waveform.riseTimeArray = riseTimeArray;
    waveform.tau_m = 10;
    filtering = 'on';
end


%% Detect events:
if isfield(detectionParameters,'voltageClamp') && detectionParameters.voltageClamp
    [minis, filterV, spectrum, waveform, F] = detectMinis(dataProperties.current, excludedTimes, detectionParameters, filtering, waveform, parallelCores, options);
    %[minis, filterV, spectrum, waveform, F, smoothV] = detectMinis(dataProperties.current, excludedTimes, detectionParameters, filtering, waveform, parallelCores, options);
    %writeABF(single([dataProperties.sweep; smoothV]), ['D:\PhD\previous\Guy_Major\Data\Voltage_clamp' filesep 'p131a_0015_sw7,8,0016_sw5,10,0017sw1_smoothed.abf'], 1000/dataProperties.dt, {'mV';'pA'});
else
    [minis, filterV, spectrum, waveform, F] = detectMinis(dataProperties.sweep, excludedTimes, detectionParameters, filtering, waveform, parallelCores, options);
    %[minis, filterV, spectrum, waveform, F, smoothV] = detectMinis(dataProperties.sweep, excludedTimes, detectionParameters, filtering, waveform, parallelCores, options);
    %writeABF(single([smoothV; dataProperties.current]), [dataProperties.filename '_smoothed.abf'], 1000/dataProperties.dt, {'mV';'pA'});
end


%% Display the event histograms:
if display
    [amplitudes1D, riseTimes1D, events2D] = classifyMinis(minis(:,4), minis(:,12), classificationParameters);
    
    % Drawing uni-dimensional fits:
    f1 = figure('position', [50 50 1200 600]);
    figure(f1);
    subplot(2,1,1,'replace');
    hold on
    targetAmpcs = cumsum(amplitudes1D)/sum(amplitudes1D);
    targetAmp99prc = find(targetAmpcs > .99, 1);
    iEndAmp = min([targetAmp99prc+1 length(amplitudeArray)]);
    xlim([amplitudeArray(1) amplitudeArray(iEndAmp)]);
    plot([0 amplitudeArray(end)], [0 0], 'k-');
    plot(amplitudeArray, amplitudes1D, 'b.-');
    set(f1, 'NumberTitle', 'off');
    set(f1, 'Name', 'One-dimensional histograms');
    title('Amplitude distribution');
    xlabel('Amplitude(mV)');
    ylabel('Number of Events');
    hold off
    
    subplot(2,1,2,'replace');
    hold on
    targetRTcs = cumsum(riseTimes1D)/sum(riseTimes1D);
    targetRT99prc = find(targetRTcs > .99, 1);
    iEndRT = min([targetRT99prc+1 length(riseTimeArray)]);
    xlim([riseTimeArray(1) riseTimeArray(iEndRT)]);
    plot([0 riseTimeArray(end)], [0 0], 'k-');
    plot(riseTimeArray, riseTimes1D,'b.-');
    if strcmpi(detectionParameters.RTinterval, '10-90%')
        title('10-90% rise time distribution');
        xlabel('10-90% rise times (ms)');
    elseif strcmpi(detectionParameters.RTinterval, '20-80%')
        title('20-80% rise time distribution');
        xlabel('20-80% rise times (ms)');
    end
    ylabel('Number of Events');
    hold off
    
    
    %% Drawing two-dimensional fits:
    f2 = figure('position', [50 50 1200 600]);
    figure(f2);
    hold on
    xlim([amplitudeArray(1)-ampLim2D amplitudeArray(iEndAmp)+ampLim2D]);
    ylim([riseTimeArray(1)-RTLim2D riseTimeArray(iEndRT)+RTLim2D]);
    imagesc(amplitudeArray, riseTimeArray, events2D);
    set(gca,'YDir','normal');
    set(f2, 'NumberTitle', 'off');
    set(f2, 'Name', 'Two-dimensional histogram');
    title('2D event distribution');
    xlabel('Amplitude(mV)');
    if strcmpi(detectionParameters.RTinterval, '10-90%')
        ylabel('10-90% rise times (ms)');
    elseif strcmpi(detectionParameters.RTinterval, '20-80%')
        ylabel('20-80% rise times (ms)');
    end
    colorbar;
    hold off
    
    F = [f1 f2 F];
end