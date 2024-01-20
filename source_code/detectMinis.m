function [minis, filterV, spectrum, waveform, varargout] = detectMinis(V, excludedT, searchParameters, filtering, varargin)
% DETECTMINIS finds miniature excitatory postsynaptic potential-like
% (mEPSPs-like or minis-like) events, their locations and other associated
% characteristics in the supplied electrophysiological recording data
% (either voltage or current clamp).
%
%   MINIS = DETECTMINIS(V, excludedT, searchParameters, filtering)
%   finds minis-like events, their locations and other associated
%   characteristics. EXCLUDEDT is the vector of excluded times.
%   SEARCHPARAMETERS is a structure variable with fields:
%       'sampleInterval' - the sampling unit size in miliseconds (ms);
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
%   FILTERING is a structure variable that is used by the filterMinis
%   function. It contains four fields:
%       'state' - a string that controls the band-stop filtering of the
%           electrophysiological data with values 'on' and 'off'
%           corresponding to filtering and leaving the data unfiltered,
%           respectively. The third option 'spectrum' returns frequency
%           spectrum vector without filtering the data.
%       'nSweeps' - the number of sweeps concatenated into a single
%           recording data file.
%       'pulseEnd' - a scalar corresponding to the end time of the last
%           pulse.
%       'filtfs' - frequencies to be filtered out (optional). It is a cell
%           array of the following example form: {'50, 150'}.
%   The output variable MINIS is a matrix with rows corresponding to the
%   detected minis-like events sorted according their temporal order and
%   columns (from left to right) corresponding to the following
%   characteristics of detected event:
%       (1)  - the peak value, mV or nA;
%       (2)  - the peak time, ms;
%       (3)  - the peak index;
%       (4)  - the amplitude, mV;
%       (5)  - the baseline value, mV;
%       (6)  - the start of a baseline, ms;
%       (7)  - the end of a baseline, ms;
%       (8)  - the element-wise start of a baseline, ms;
%       (9)  - the element-wise end of a baseline, ms;
%       (10) - the element-wise 0-100% rise time (from the start of the
%              baseline to the peak).
%       (11) - the 0-100% rise time, ms.
%       (12) - the 10-90% rise time, ms.
%       (13) - the 10% rise time mark, ms.
%       (14) - the 50% rise time mark, ms.
%       (15) - the 90% rise time mark, ms.
%       (16) - the index of the 10% rise time mark.
%       (17) - the index of the 50% rise time mark.
%       (18) - the index of the 90% rise time mark.
%       (19) - membrane potential or current value at the 10% rise time
%              mark, mV or nA.
%       (20) - membrane potential or current value at the 50% rise time
%              mark, mV or nA.
%       (21) - membrane potential or current value at the 90% rise time
%              mark, mV or nA.
%       (22) - 1/e decay, ms.
%
%   [MINIS, FILTERV] = DETECTMINIS(...)
%   In addition returns a band-stop/notch filtered electrophysiological
%   recording vector FILTERV. This option is useful when minis-like event
%   detection is performed as part of the manual or automatic distribution
%   fitting using the 'minis' program. The output vector FILTERV is already
%   filtered and so can be later re-used without filtering.
%
%   [MINIS, FILTERV, SPECTRUM] = DETECTMINIS(...)
%   Returns the amplitude, power, and phase spectra of the recording after
%   a fast Fourier transform if the FILTERING state is set to 'on' or
%   'spectrum'. Otherwise returns an empty array. SPECTRUM is a matrix
%   composed of four row vectors. The first vector contains amplitudes. The
%   second vector contains a power spectrum which should correspond to
%   decibels (dB). The third vector contains phase information. The fourth
%   vector contains frequencies.
%
%   [MINIS, FILTERV, SPECTRUM, WAVEFORM] = DETECTMINIS(..., waveform)
%   WAVEFORM is a structure variable containing the following fields:
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
%   The output variable WAVEFORM is a structure containing the fields:
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
%   If the input waveform variable is not set to 1 or 2, the output
%   waveform variable is empty.
%
%   MINIS = DETECTMINIS(..., waveform, parallelCores)
%   performs the detection algorithm in a parallel fashion. PARALLELCORES
%   is a scalar variable corresponding to the number of processor cores
%   allocated for the task. The number of available cores depends on the
%   your hardware. In order to run the function in parallel, at least two
%   Matlab workers/labs have to be available. For opening additional matlab
%   workers see the documentation for the "matlabpool" function (or type in
%   the command line "help matlabpool").
%
%   MINIS = DETECTMINIS(..., waveform, parallelCores, options)
%   in addition allows a user to control the visualisation of the minis
%   detection process. OPTIONS is a structure variable containing the
%   following fields:
%       'interactMode' - either true or false. Choose true for visualising
%           the detection process and receiving prompts for data inspection.
%           Otherwise choose false. False is the default option.
%           Interactive mode is unavailable in parallel processing due to
%           the restrictions of the 'parfor' loop.
%       'pauseTime' - the minimum waiting time in the interactive mode
%           given for performance inspection.
%       'filename' - the name of the data file.
%       'dataType' - the data type (i.e., 'Membrane potential' or 'Current');
%       'dataUnits' - data measurement units (i.e., 'mV' or 'nA');
%       'summaryPlot' - either true or false. Choose true to obtain a
%           summary plot with detection results displayed, as well as the
%           Fourier time-to-frequency transform; or choose false otherwise.
%           The default state is false.
%       'edit' - a logical field that controls whether to allow the manual
%           editing of detected events. True is to allow and false is not.
%           The deffault value is false.
%       'smooth' - a logical field with true representing smoothing of the
%           electrophysiological recording data and false otherwise. The
%           default value is true.
%       'SD' - a vector variable controlling the estimation of the standard
%           deviation (SD) of the recording trace. The vector contains
%           four elements each of which could be set to either 1 (true) or
%           0 (false). The following combinations are used to calculate
%           different types of SD:
%               [1 0 0 0] - calculates the SD of the entire recording trace
%                   without any adjustments for non-stationarity of the
%                   data.
%               [0 1 0 0] - calculates the SD of the entire recording trace
%                   with non-stationarity minimised. The calculation is
%                   based on 15ms-window average.
%               [0 0 1 0] - calculates the SD of the entire recording trace
%                   with non-stationarity minimised. The calculation is
%                   based on 30ms-window average.
%               [0 0 0 1] - calculates the SD of the intervals containing
%                   no detected minis-like events. This calculation
%                   estimates the SD of the pseudo noise and is based on
%                   15ms-window average.
%               Other combinations are possible allowing to calculate more
%                   than one type of SD at the same time.
%           The default combination is [0 0 0 0] where no SD calculations
%           are carried out. If the user specifies any of the combinations
%           containing ones, the estimates are appended to the MINIS output
%           array as additional collumns 22-25 in the order they are listed.
%
%   [MINIS, FILTERV, SPECTRUM, WAVEFORM, F] = DETECTMINIS(...)
%   also returns a handle F to a figure with the recording trace and
%   detected minis (the summary plot), as well as, frequency power spectrum
%   figure after Fourier transform (if filtering was an option).
%
%   [MINIS, FILTERV, SPECTRUM, WAVEFORM, F, SMOOTHV] = DETECTMINIS(...)
%   In addition returns a smoothed electrophysiological recording vector
%   SMOOTHV.
%



%% Check if program termination was requested
global STOP
if STOP
  disp('Program execution terminated');
  error('Interrupted by user');
end


%% Identify input variables:
if nargin >= 5
    waveform = varargin{1};
    if nargin >= 6
        parallelCores = varargin{2};
        parallelCores = round(parallelCores);
        if parallelCores <= 0
            parallelCores = 1;
        end
    end
else
    waveform.estimate = 0;
    parallelCores = 1;
end

if nargin == 7
    options = varargin{3};
    if isfield(options,'interactMode')
        interactMode = options.interactMode;
        if interactMode
            if isfield(options,'pauseTime')
                pauseTime = options.pauseTime;
            else
                pauseTime = .5;
            end
        end
    else
        interactMode = false;
    end
    if isfield(options,'filename')
        filename = options.filename;
        if isstruct(filename) && isfield(filename, 'filename')
            filename = filename.filename;
        end
    else
        filename = 'a data file';
    end
    if isfield(options,'dataType') && (strcmpi(options.dataType,'Membrane potential') || strcmpi(options.dataType,'Current'))
        dataType = options.dataType;
        if strcmpi(dataType,'Membrane potential')
            dataUnits = 'mV'; %#ok<NASGU>
        else
            dataUnits = 'nA'; %#ok<NASGU>
        end
    else
        dataType = 'Membrane potential or current';
    end
    if isfield(options,'dataUnits') && (strcmpi(options.dataUnits,'mV') || strcmpi(options.dataUnits,'nA'))
        dataUnits = options.dataUnits;
    else
        dataUnits = '(mV or nA)';
    end
    if isfield(options,'summaryPlot') && options.summaryPlot
        summaryPlot = true;
    else
        summaryPlot = false;
    end
    if isfield(options,'edit') && options.edit
        edit = true;
    end
    if (isfield(options,'smooth') && ~options.smooth) ||...
            (isfield(searchParameters,'smoothWindow') && ~searchParameters.smoothWindow)
        smooth = false;
    else
        smooth = true;
    end
    if isfield(options,'SD')
        SD = logical(options.SD);
    else
        SD = [false false false false];
    end
else
    interactMode = false;
    summaryPlot = false;
    edit = false;
    smooth = true;
    SD = [false false false false];
end
clear options varargin

dt = searchParameters.sampleInterval;                                       % The sample interval, ms
initt = dt*(0:length(V) - 1);                                               % The time vector, ms
smoothWindows = [searchParameters.smoothWindow searchParameters.smoothWindowLite];



%% Check unique device signature
if ~fafTest
    if ispc
      [~,uniqueDeviceSignature] = system('wmic diskdrive get serialnumber');
    elseif ismac
      [~,uniqueDeviceSignature] = system('system_profiler SPHardwareDataType | grep Serial');
    elseif isunix
      [~,uniqueDeviceSignature] = system('ls -la /dev/disk/by-uuid');
    end
    %if ~contains(uniqueDeviceSignature, '2G2720') || ~contains(uniqueDeviceSignature, '064419') % MD1
    %if ~contains(uniqueDeviceSignature, 'ACE4_2E00_2520_4A01_') || ~contains(uniqueDeviceSignature, '2EE4_AC00_0000_0001') % MD2
    %if ~contains(uniqueDeviceSignature, 'ACE4_2E00_15FF_E770_') || ~contains(uniqueDeviceSignature, '2EE4_AC00_0000_0001') % GM
    %if ~contains(uniqueDeviceSignature, 'E823_8FA6_BF53_0001_') || ~contains(uniqueDeviceSignature, '001B_448B_497C_2616') % UID
    if false
        minis = []; filterV = []; spectrum = struct();
        waveform.parameters.averageAmp = [];
        waveform.parameters.medianAmp = [];
        waveform.parameters.tau_m = []; waveform.estimate = [];
        waveform.tau_m = []; waveform.nSweeps = [];
        waveform.riseTimeArray = []; varargout{1} = [];
        return
    end
end



%% Invert data
if isfield(searchParameters,'downGoing') && searchParameters.downGoing
    V = -V;
end



%% Band-stop filter the regular (systemic) noise:
if strcmpi(filtering.state,'on')
    if summaryPlot
        if exist('filename','var')
            [V, spectrum, f1] = filterMinis(V, dt, filtering, summaryPlot, filename);
        else
            [V, spectrum, f1] = filterMinis(V, dt, filtering, summaryPlot);
        end
    else
        [V, spectrum] = filterMinis(V, dt, filtering, summaryPlot);
    end
elseif strcmpi(filtering.state,'spectrum')
    [~, spectrum] = filterMinis(V, dt, filtering, summaryPlot);
else
    spectrum = [];
end
filterV = V;



%% Detection Algorithm:
% Detect initial peaks:
if smooth
    initvsmooth = gsmooth(V, smoothWindows(1), 3);
else
    initvsmooth = V;
end
[~, origPeakIndexArray] = findpeaks(double(initvsmooth));
clear smooth filter

% Calculate excluded times:
excludedT(excludedT < 0 & excludedT > initt(end)) = [];
excludedIndices = round(excludedT/dt) + 1;
exclIndLogical = zeros(1, length(initt));
exclIndLogical(excludedIndices) = 1;
circ = searchParameters.SWstart + round(searchParameters.BLduration/dt);
response = ones(1, 1 + 2*circ);
excludedIndices = logical(conv(exclIndLogical, response, 'same'));
excludedIndices(1 : (searchParameters.SWstart + searchParameters.refractoryPeriod)/dt - 1) = true;
excludedT = initt(excludedIndices);

% Remove peaks falling within excluded times:
peaksToKeep = zeros(1, length(initt));
peaksToKeep(origPeakIndexArray) = 1;
peaksToKeep(logical(peaksToKeep) & excludedIndices) = 0;
initPeakIndexArray = find(logical(peaksToKeep));

% Plot membrane potential data:
if interactMode
    inspectSmoothing = true;                                                % Set to true if you want to inspect the smoothing of data
    if inspectSmoothing
        vplotf2 = [V; initvsmooth];
        nameStringf2 = 'Examine smoothing: Raw vs. smoothed data';
        titleStringf2 = sprintf('Data from %s', filename);
        optionsf2 = struct('nameString', nameStringf2, 'titleString', titleStringf2, 'dataType', dataType, 'dataUnits', dataUnits);
        f2 = plotData(initt, vplotf2, optionsf2);
        
        vsmoothLight = gsmooth(V, smoothWindows(2), 3);
        vplotf3 = [vsmoothLight; initvsmooth];
        nameStringf3 = 'Examine smoothing: Lightly vs. normally smoothed data';
        titleStringf3 = sprintf('Smoothed data from %s', filename);
        optionsf3 = struct('nameString', nameStringf3, 'titleString', titleStringf3, 'dataType', dataType, 'dataUnits', dataUnits);
        f3 = plotData(initt, vplotf3, optionsf3);
        pause
        if ishandle(f2)
            close(f2)
        end
        if ishandle(f3)
            close(f3)
        end
    end
    
    % Initialise the figure for visualising minis detection:
    nameStringf4 = 'Examine events: Smoothed data only';
    titleStringf4 = sprintf('Smoothed data from %s', filename);
    optionsf4 = struct('nameString', nameStringf4, 'titleString', titleStringf4, 'dataType', dataType, 'dataUnits', dataUnits);
    f4 = plotData(initt, initvsmooth, optionsf4);
    figure(f4);
    axes3 = get(gcf,'CurrentAxes');
    hold on
    % Plot excluded times:
    vexcluded = initvsmooth(excludedIndices);
    plot(excludedT, vexcluded, 'y.', 'markersize', 5);
    pause(pauseTime*2);
end
clear V smoothWindows;

% PLOT1 - Examine the initial peaks and then dim them:
if interactMode
    peakTimes = initPeakIndexArray*dt - dt;
    initialPeaks = initvsmooth(initPeakIndexArray);
    figure(f4);
    hold on;
    p1 = plot(peakTimes, initialPeaks, '^', 'markersize', 10, 'markerfacecolor', 'r', 'markeredgecolor', 'r');
    pause
    delete(p1);
    plot(peakTimes, initialPeaks, '^', 'color', 'r');
    hold off;
end

% Detect minis-like events:
switch parallelCores
    case 1                                                                  % Serial detection
        if interactMode
            finalMinis = detectionAlgorithmSerial(initt, initvsmooth, searchParameters, initPeakIndexArray, origPeakIndexArray, interactMode, pauseTime,...
                f4, axes3);
        else
            finalMinis = detectionAlgorithmSerial(initt, initvsmooth, searchParameters, initPeakIndexArray, origPeakIndexArray);
        end
    otherwise                                                               % Parallel detection
        % Divide the recording sweep for parallel processing:
        add = 20;
        addOrig = ceil((length(origPeakIndexArray) - length(initPeakIndexArray))/parallelCores);
        initLength = ceil(length(initPeakIndexArray)/parallelCores);
        origLength = ceil(length(origPeakIndexArray)/parallelCores);
        peakIndexArray = zeros(parallelCores, initLength+add);
        originalPeakIndexArray = zeros(parallelCores, origLength+add+2*addOrig);
        peakIndexArray(1,:) = initPeakIndexArray(1:initLength+add);
        originalPeakIndexArray(1,:) = origPeakIndexArray(1:origLength+add+2*addOrig);
        peakIndexArray(end,:) = initPeakIndexArray(end-initLength-add+1:end);
        originalPeakIndexArray(end,:) = origPeakIndexArray(end-origLength-add-2*addOrig+1:end);
        for iCores = 2:parallelCores-1
            peakIndexArray(iCores,:) = initPeakIndexArray((iCores-1)*initLength-add/2+1:iCores*initLength+add/2);
            originalPeakIndexArray(iCores,:) = origPeakIndexArray((iCores-1)*origLength-add/2-addOrig+1:iCores*origLength+add/2+addOrig);
        end
        % Carry out the detection:
        finalMinis = zeros(initLength+add,22,parallelCores);
        parfor iCore = 1:parallelCores
            coreMinis = detectionAlgorithmParallel(initt, initvsmooth, searchParameters, peakIndexArray(iCore,:), originalPeakIndexArray(iCore,:));
            finalMinis(:,:,iCore) = coreMinis;
        end
end

% Concatenate the detection output:
vdim = size(finalMinis,1);
hdim = size(finalMinis,2);
tempMinis = zeros(parallelCores*vdim,hdim);
for iCore = 1:parallelCores
    tempMinis((iCore-1)*vdim+1:iCore*vdim,:) = finalMinis(:,:,iCore);
end
uniqueMinis = unique(tempMinis,'rows','first');
minis = sortrows(uniqueMinis,2);
minis(1,:) = [];
minis = minis(~minis(:,3)==0,:);
assert(~sum(minis(:,3)==0));
clear parallelCores interactMode peakIndexArray pauseTime f4 axes3 finalMinis iCore coreMinis vdim hdim tempMinis uniqueMinis

% PLOT6 - Mark detected minis for examination:
if summaryPlot || (exist('edit','var') && edit)
    nameStringf5 = 'Summary: Detected mini-like events';
    titleStringf5 = sprintf('Smoothed data from %s', filename);
    optionsf5 = struct('nameString', nameStringf5, 'titleString', titleStringf5, 'dataType', dataType, 'dataUnits', dataUnits);
    f5 = plotData(initt, initvsmooth, optionsf5);
    figure(f5);
    hold on;
    p10 = plot(minis(:,2), minis(:,1), '^', 'markersize', 10, 'markeredgecolor', 'r');
    p11 = plot(minis(:,15), minis(:,21), 'r.', 'markersize', 10);
    p12 = plot(minis(:,14), minis(:,20), 'r.', 'markersize', 10);
    p13 = plot(minis(:,13), minis(:,19), 'r.', 'markersize', 10);
    p14 = plot(minis(:,6), minis(:,5), 'c.', 'markersize', 10);
    p15 = plot(minis(:,7), minis(:,5), 'b.', 'markersize', 10);
    % Plot excluded times:
    excludedIndices = round(excludedT/dt) + 1;
    vexcluded = initvsmooth(excludedIndices);
    p16 = plot(excludedT, vexcluded, 'y.', 'markersize', 5);
    legend([p10, p11, p14, p15, p16],{'Peak','RT markers','BL start','BL end','Excluded'}, 'Location','NorthEast');
    hold off;
    
    % Manual adjustment:
    if exist('edit','var') && edit
        button = questdlg('Edit the detected events?','Edit Events','Yes','No','No');
        if strcmpi(button, 'Yes')
            bringHelp;
            axes5 = get(gcf,'CurrentAxes');
            [minis, f5] = manualAdjust(minis, initt, initvsmooth, f5, [p10, p11, p12, p13, p14, p15], axes5,...
                searchParameters.BLduration, searchParameters.RTinterval);
        end
    end
    if exist('f1','var')
        varargout{1} = [f5 f1];
    else
        varargout{1} = f5;
    end
else
    varargout{1} = 0;
end
varargout{2} = initvsmooth;
clear summaryPlot edit nameStringf5 titleStringf5 optionsf5 f1 f5 initt p10 p11 p12 p13 p14 p15 vexcluded excludedIndices button axes5 dataType dataUnits filename



%% Estimate the average waveform:
if waveform.estimate || SD(4)
    if waveform.estimate == 1
        button = questdlg('Do you want to estimate the average waveform?','Obtain Average Waveform','Yes','No','Yes');
        if strcmpi(button, 'Yes')
            [waveform.averageTrace, averageAmp, medianAmp, waveform.riseTimeDist, waveform.parameters, waveform.F] = averageWaveEffect(minis,...
                initvsmooth, dt, waveform.tau_m, excludedT, 10, searchParameters.RTinterval, waveform.riseTimeArrayExt, true);
           waveform.parameters.averageAmp = averageAmp;
           waveform.parameters.medianAmp = medianAmp;
        else
            waveform = [];
        end
    elseif waveform.estimate == 3
        [waveform.averageTrace, averageAmp, medianAmp, waveform.riseTimeDist, waveform.parameters, waveform.F] = averageWaveEffect(minis, initvsmooth, dt, waveform.tau_m,...
                excludedT, 10, searchParameters.RTinterval, waveform.riseTimeArrayExt, true);
        waveform.parameters.averageAmp = averageAmp;
        waveform.parameters.medianAmp = medianAmp;
    else
        [waveform.averageTrace, averageAmp, medianAmp, waveform.riseTimeDist, waveform.parameters] = averageWaveEffect(minis, initvsmooth, dt, waveform.tau_m, excludedT,...
            10, searchParameters.RTinterval, waveform.riseTimeArrayExt, false);
        waveform.parameters.averageAmp = averageAmp;
        waveform.parameters.medianAmp = medianAmp;
    end
else
    waveform = [];
end



%% Estimate the standard deviation:
if any(SD)
    STD = [0 0 0 0];
    
    % Total non-stationarity unadjusted SD:
    if SD(1)
        STD(1) = std(initvsmooth(~exclIndLogical),1);
    end
    
    if SD(2) || SD(4)
        window = round(15/dt);
        rows = floor(length(initvsmooth)/window);
        vert = repmat(window*(0:rows-1)',1,window);
        horz = repmat((1:window),rows,1);
        inds = vert + horz;
        arrayV = initvsmooth(inds);
        
        % Non-stationarity adjusted SD (15ms window):
        if SD(2)
            arrayExcl = logical(sum(exclIndLogical(inds),2));
            totSTD = std(arrayV,1,2);
            STD(2) = mean(totSTD(~arrayExcl));
        end
        
        % Estimated pseudo noise SD:
        if SD(4)
            peakExclusion = zeros(length(initvsmooth), 1);
            peakExclusion(minis(:,3)) = 1;
            if isempty(waveform) || isempty(waveform.parameters.tau_m)
                tau_m = 10;
            else
                tau_m = waveform.parameters.tau_m;
            end
            response = ones(1, round((6*tau_m)/dt));
            peakExclusion = logical(conv(peakExclusion', response, 'same'));
            exclIndLogicalNoise = logical(logical(exclIndLogical) + peakExclusion);
            arrayExcl = logical(sum(exclIndLogicalNoise(inds),2));
            totSTD = std(arrayV,1,2);
            STD(4) = mean(totSTD(~arrayExcl));
        end
    end
    
    % Non-stationarity adjusted SD (30ms window):
    if SD(3)
        window = round(30/dt);
        rows = floor(length(initvsmooth)/window);
        vert = repmat(window*(0:rows-1)',1,window);
        horz = repmat((1:window),rows,1);
        inds = vert + horz;
        arrayV = initvsmooth(inds);
        
        arrayExcl = logical(sum(exclIndLogical(inds),2));
        totSTD = std(arrayV,1,2);
        STD(3) = mean(totSTD(~arrayExcl));
    end
    
    if ~isempty(STD)
        minis = [minis repmat(STD(STD ~= 0), size(minis,1), 1)];
    end
end
end










%% Local functions:
function coreMinis = detectionAlgorithmSerial(t, vsmooth, searchParameters, peakIndices, originalPeakInd, varargin)
% DETECTIONALGORITHMSERIAL is a serial minis detection algorithm and is a
% helper subfunction of DETECTMINIS.
%



% Assign major variables:
if nargin > 5
    interactMode = varargin{1};
    pauseTime = varargin{2};
    f4 = varargin{3};
    axes3 = varargin{4};
else
    interactMode = false;
end
RTtype = searchParameters.RTinterval;                                       % The rise time type (20-80% or 10-90%)
dt = searchParameters.sampleInterval;                                       % The sample interval, ms
SWupbound = searchParameters.SWstart;                                       % The beginning of the search window, ms
nSWupbound = floor(SWupbound/dt);                                           % The element-wise beginning of the search window
tBLduration = searchParameters.BLduration;                                  % The baseline duration, ms
nBLduration = round(tBLduration/dt);                                        % The element-wise baseline duration
tRefractory = searchParameters.refractoryPeriod;                            % The refractory period, ms
nRefractory = ceil(tRefractory/dt);                                         % The element-wise refractory period
Amplobound = searchParameters.Amplobound;                                   % The lower bound on the amplitude size, mV
Ampupbound = searchParameters.Ampupbound;                                   % The upper bound on the amplitude size, mV
lengthTempStore = length(peakIndices);


% Initialise storage variables:
peakCounter = 0;                                                            % Established peak counter variable
finalisedPeaks = zeros(lengthTempStore,1);                                  % The array of confirmed peaks
peakTimes = zeros(lengthTempStore,1);
finalisedPeakIndices = zeros(lengthTempStore,1);                            % The indices of confirmed peaks
amplitudes = zeros(lengthTempStore,1);                                      % and their associated characteristics (below)
baselines = zeros(lengthTempStore,1);
tBaselines = zeros(lengthTempStore,2);
nBaselines = zeros(lengthTempStore,2);
elementRiseTimes = zeros(lengthTempStore,1);
riseTimes = zeros(lengthTempStore,1);
riseTimes1090 = zeros(lengthTempStore,1);
t10s = zeros(lengthTempStore,1);
t50s = zeros(lengthTempStore,1);
t90s = zeros(lengthTempStore,1);
n10s = zeros(lengthTempStore,1);
n50s = zeros(lengthTempStore,1);
n90s = zeros(lengthTempStore,1);
v10s = zeros(lengthTempStore,1);
v50s = zeros(lengthTempStore,1);
v90s = zeros(lengthTempStore,1);
decays = NaN(lengthTempStore,1);


% Major loop:
for iPeak = 1:length(peakIndices)
    peakIndex = peakIndices(iPeak);
    peakTime = t(peakIndex);
    
    % Filter out peaks that do not meat the inclusion criteria:
    % 1. Test whether the local maximum (peak) is the highest peak in the vicinity ahead:
    iEndRefractory = peakIndex + nRefractory;
    refractoryPeakIndices = peakIndices(peakIndices > peakIndex & peakIndices <= iEndRefractory);
    refractoryPeaks = vsmooth(refractoryPeakIndices);
    maxRefractoryPeak = max(refractoryPeaks);
    if isempty(maxRefractoryPeak) || vsmooth(peakIndex) > maxRefractoryPeak
        
        % 2. Test whether the local maximum (peak) has no established peaks in the vicinity behind:
        if peakCounter && peakTime - peakTimes(peakCounter) > tRefractory || ~peakCounter
            
            % PLOT2.1 - Mark a potential mini:
            if interactMode
                figure(f4);
                hold on;
                set(axes3,'XLim',[t(peakIndex - nSWupbound) t(peakIndex + nRefractory)]);
                p2 = plot(peakTime, vsmooth(peakIndex), '^', 'markersize', 10, 'markerfacecolor', 'r', 'markeredgecolor', 'r');
                hold off;
                pause(pauseTime);
            end
            
            % 3. Test whether the local maximum (peak) is in the interval of acceptable amplitudes:
            % 3.1. Find the baseline of the peak:
            if peakCounter
                iSWstart = max([peakIndex - nSWupbound + 1, finalisedPeakIndices(peakCounter)]);
            else
                iSWstart = max([peakIndex - nSWupbound + 1, round(.5*nSWupbound)]);
            end
            [trough, iTrough] = min(vsmooth(iSWstart:peakIndex));
            iTrialBLend = min([iSWstart - 1 + iTrough + round(.2*nBLduration) peakIndex-4]);
            if peakCounter
                iTrialBLstart = max([iTrialBLend - nBLduration + 1 finalisedPeakIndices(peakCounter)]);
            else
                iTrialBLstart = iTrialBLend - nBLduration + 1;
            end
            trialBL = mean(vsmooth(iTrialBLstart:iTrialBLend));
            
            % PLOT3.1 - Bracket out the baseline:
            if interactMode
                figure(f4);
                hold on;
                p3 = plot(t(iTrialBLstart), trialBL, 'c.', 'markersize', 10);
                p4 = plot(t(iTrialBLend), trialBL, 'b.', 'markersize', 10);
                hold off;
                pause(.5*pauseTime);
            end
            
            % 3.2. Test whether there are no higher peaks between the baseline and the proper peak (for removing edge effects):
            if isempty(vsmooth(originalPeakInd >= iSWstart - 1 + iTrough & originalPeakInd < peakIndex - 1))...
                    || ~(max(vsmooth(originalPeakInd(originalPeakInd >= iSWstart - 1 + iTrough & originalPeakInd < peakIndex - 1))) >= vsmooth(peakIndex))
                
                % 3.3. Test the amplitude:
                amp = vsmooth(peakIndex) - trialBL;
                if amp > Amplobound && amp < Ampupbound
                    
                    % 4. Test the rise time spread:
                    if interactMode
                        [t1090, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, t(iTrialBLend:peakIndex), vsmooth(iTrialBLend:peakIndex), trialBL,...
                            interactMode, pauseTime);
                    else
                        [t1090, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, t(iTrialBLend:peakIndex), vsmooth(iTrialBLend:peakIndex), trialBL);
                    end
                    part1RT = i50 - i10;
                    part2RT = i90 - i50;
                    
                    if part1RT/part2RT >= 5                                 % In case the first half is considerably larger than the second one
                        % 4.1. Correct the baseline if possible and re-test the amplitude:
                        newSWstart = peakIndices(find(peakIndices >= iTrialBLend, 1));
                        if newSWstart < peakIndex
                            if interactMode
                                [pass, amp, trialBL, iTrialBLstart, iTrialBLend, t1090, t10, i10, t50, i50, t90, i90, p3, p4] = correctBL(vsmooth, t,...
                                    RTtype, newSWstart, peakIndex, nBLduration, Amplobound, Ampupbound, peakCounter, finalisedPeakIndices, interactMode,...
                                    pauseTime, f4, p3, p4);
                                p3 = cell2mat(p3);
                                p4 = cell2mat(p4);
                            else
                                [pass, amp, trialBL, iTrialBLstart, iTrialBLend, t1090, t10, i10, t50, i50, t90, i90] = correctBL(vsmooth, t, RTtype,...
                                    newSWstart, peakIndex, nBLduration, Amplobound, Ampupbound, peakCounter, finalisedPeakIndices, interactMode);
                            end
                            if pass
                                part1RT = i50 - i10;
                                part2RT = i90 - i50;
                                
                                if part1RT/part2RT >= 5                     % In case the first half is still considerably larger than the second one
                                    % 4.2. Repeat correction:
                                    newSWstart = peakIndices(find(peakIndices >= iTrialBLend, 1));
                                    if newSWstart ~= peakIndex
                                        if interactMode
                                            [pass, amp, trialBL, iTrialBLstart, iTrialBLend, t1090, t10, i10, t50, i50, t90, i90, p3, p4] = correctBL(...
                                                vsmooth, t, RTtype, newSWstart, peakIndex, nBLduration, Amplobound, Ampupbound, peakCounter,...
                                                finalisedPeakIndices, interactMode, pauseTime, f4, p3, p4);
                                            p3 = cell2mat(p3);
                                            p4 = cell2mat(p4);
                                        else
                                            [pass, amp, trialBL, iTrialBLstart, iTrialBLend, t1090, t10, i10, t50, i50, t90, i90] = correctBL(vsmooth, t,...
                                                RTtype, newSWstart, peakIndex, nBLduration, Amplobound, Ampupbound, peakCounter, finalisedPeakIndices,...
                                                interactMode);
                                        end
                                    end
                                    if ~pass
                                        % PLOT5 - Destroy unsuccessful mini candidate:
                                        if interactMode
                                            figure(f4);
                                            hold on;
                                            delete([p3, p4]);
                                            p5 = plot(peakTime, vsmooth(peakIndex), 'x', 'markersize', 20, 'color', 'k', 'LineWidth', 2);
                                            pause(2*pauseTime);
                                            delete(p2);
                                            delete(p5);
                                            hold off;
                                        end
                                        continue
                                    end
                                end
                                
                            else
                                % PLOT5 - Destroy unsuccessful mini candidate:
                                if interactMode
                                    figure(f4);
                                    hold on;
                                    delete([p3, p4]);
                                    p5 = plot(peakTime, vsmooth(peakIndex), 'x', 'markersize', 20, 'color', 'k', 'LineWidth', 2);
                                    pause(2*pauseTime);
                                    delete(p2);
                                    delete(p5);
                                    hold off;
                                end
                                continue
                            end
                        end
                        
                    elseif part1RT/part2RT <= .2                            % In case the first half is considerably smaller than the second one
                        % 4.3. Shift the peak backwards if possible:
                        newPeakIndex = peakIndices(find(peakIndices < peakIndex, 1, 'last'));
                        if ~isempty(newPeakIndex) && newPeakIndex > iTrialBLend
                            newPeakTime = t(newPeakIndex);
                            
                            % PLOT2.2 - Remark the peak:
                            if interactMode
                                figure(f4);
                                hold on;
                                delete(p2);
                                set(axes3,'XLim',[t(newPeakIndex - nSWupbound) t(newPeakIndex + nRefractory)]);
                                p2 = plot(newPeakTime, vsmooth(newPeakIndex), '^', 'markersize', 10, 'markerfacecolor', 'r', 'markeredgecolor', 'r');
                                hold off;
                                pause(pauseTime);
                            end
                            
                            % 4.4. Re-test the amplitude:
                            newAmp = vsmooth(newPeakIndex) - trialBL;
                            if newAmp > Amplobound && newAmp < Ampupbound
                                if interactMode
                                    [t1090, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, t(iTrialBLend:newPeakIndex),...
                                        vsmooth(iTrialBLend:newPeakIndex), trialBL, interactMode, pauseTime);
                                else
                                    [t1090, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, t(iTrialBLend:newPeakIndex),...
                                        vsmooth(iTrialBLend:newPeakIndex), trialBL);
                                end
                                peakIndex = newPeakIndex;
                                peakTime = newPeakTime;
                                amp = newAmp;
                            else
                                % PLOT2.3 - Remark the peak again:
                                if interactMode
                                    figure(f4);
                                    hold on;
                                    delete(p2);
                                    set(axes3,'XLim',[t(peakIndex - nSWupbound) t(peakIndex + nRefractory)]);
                                    p2 = plot(peakTime, vsmooth(peakIndex), '^', 'markersize', 10, 'markerfacecolor', 'r', 'markeredgecolor', 'r');
                                    hold off;
                                    pause(pauseTime);
                                end
                            end
                        end
                    end
                    
                    % 5. Check whether the baseline does not deviate too far from the 10% rise time mark:
                    if (strcmpi(RTtype, '10-90%') && iTrialBLstart - 1 + i10 - iTrialBLend > ceil(.5*nBLduration))...
                            || (strcmpi(RTtype, '20-80%') && iTrialBLstart - 1 + i10 - iTrialBLend > nBLduration)
                        
                        % 5.1. Correct the baseline:
                        iSWstart = iTrialBLstart - 1 + i10 - ceil(.5*nBLduration);
                        [~, iTrough] = min(vsmooth(iSWstart:peakIndex));
                        iTrialBLend = min([iSWstart - 1 + iTrough + round(.2*nBLduration) peakIndex-4]);
                        iTrialBLstart = iTrialBLend - nBLduration + 1;
                        trialBL = mean(vsmooth(iTrialBLstart:iTrialBLend));
                        
                        % 5.2. Test the new amplitude:
                        amp = vsmooth(peakIndex) - trialBL;
                        if amp > Amplobound && amp < Ampupbound
                            if interactMode
                                [t1090, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, t(iTrialBLend:peakIndex), vsmooth(iTrialBLend:peakIndex), trialBL,...
                                    interactMode, pauseTime);
                            else
                                [t1090, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, t(iTrialBLend:peakIndex), vsmooth(iTrialBLend:peakIndex), trialBL);
                            end
                            
                            % PLOT3.3 - Correct the baseline:
                            if interactMode
                                figure(f4);
                                hold on;
                                delete(p3,p4);
                                pause(.5*pauseTime);
                                p3 = plot(t(iTrialBLstart), trialBL, 'c.', 'markersize', 10);
                                p4 = plot(t(iTrialBLend), trialBL, 'b.', 'markersize', 10);
                                hold off;
                                pause(.5*pauseTime);
                            end
                            
                            peakCounter = peakCounter + 1;                  % A peak has been established
                            finalisedPeaks(peakCounter) = vsmooth(peakIndex);
                            peakTimes(peakCounter) = peakTime;
                            finalisedPeakIndices(peakCounter) = peakIndex;
                            amplitudes(peakCounter) = amp;
                            baselines(peakCounter) = trialBL;
                            tBaselines(peakCounter,:) = [t(iTrialBLstart), t(iTrialBLend)];
                            nBaselines(peakCounter,:) = [iTrialBLstart, iTrialBLend];
                            elementRiseTimes(peakCounter) = peakIndex - iTrialBLstart + 1;
                            riseTimes(peakCounter) = elementRiseTimes(peakCounter)*dt;
                            riseTimes1090(peakCounter) = t1090;
                            t10s(peakCounter) = t10;
                            t50s(peakCounter) = t50;
                            t90s(peakCounter) = t90;
                            n10s(peakCounter) = iTrialBLend - 1 + i10;
                            n50s(peakCounter) = iTrialBLend - 1 + i50;
                            n90s(peakCounter) = iTrialBLend - 1 + i90;
                            v10s(peakCounter) = vsmooth(n10s(peakCounter));
                            v50s(peakCounter) = vsmooth(n50s(peakCounter));
                            v90s(peakCounter) = vsmooth(n90s(peakCounter));
                            decays(peakCounter) = calcDecay(vsmooth, baselines(peakCounter), peakIndex, dt);
                            
                            % PLOT4 - Indicate that the peak has been just established:
                            if interactMode
                                figure(f4);
                                hold on;
                                plot(t90s(peakCounter), vsmooth(n90s(peakCounter)), 'r.', 'markersize', 10);
                                plot(t50s(peakCounter), vsmooth(n50s(peakCounter)), 'r.', 'markersize', 10);
                                plot(t10s(peakCounter), vsmooth(n10s(peakCounter)), 'r.', 'markersize', 10);
                                for i = 1:3
                                    p5 = plot(peakTime, vsmooth(peakIndex), '^', 'markersize', 20, 'markerfacecolor', 'r', 'markeredgecolor', 'r');
                                    pause(pauseTime);
                                    delete(p5);
                                    pause(pauseTime);
                                end
                                hold off;
                            end
                        end
                        % Coda: If the local maximum fails in any of the 5 tests, then it is rejected.
                        
                    else
                        peakCounter = peakCounter + 1;                      % A peak has been established
                        finalisedPeaks(peakCounter) = vsmooth(peakIndex);
                        peakTimes(peakCounter) = peakTime;
                        finalisedPeakIndices(peakCounter) = peakIndex;
                        amplitudes(peakCounter) = amp;
                        baselines(peakCounter) = trialBL;
                        tBaselines(peakCounter,:) = [t(iTrialBLstart), t(iTrialBLend)];
                        nBaselines(peakCounter,:) = [iTrialBLstart, iTrialBLend];
                        elementRiseTimes(peakCounter) = peakIndex - iTrialBLstart + 1;
                        riseTimes(peakCounter) = elementRiseTimes(peakCounter)*dt;
                        riseTimes1090(peakCounter) = t1090;
                        t10s(peakCounter) = t10;
                        t50s(peakCounter) = t50;
                        t90s(peakCounter) = t90;
                        n10s(peakCounter) = iTrialBLend - 1 + i10;
                        n50s(peakCounter) = iTrialBLend - 1 + i50;
                        n90s(peakCounter) = iTrialBLend - 1 + i90;
                        v10s(peakCounter) = vsmooth(n10s(peakCounter));
                        v50s(peakCounter) = vsmooth(n50s(peakCounter));
                        v90s(peakCounter) = vsmooth(n90s(peakCounter));
                        decays(peakCounter) = calcDecay(vsmooth, baselines(peakCounter), peakIndex, dt);
                        
                        % PLOT4 - Indicate that the peak has been just established:
                        if interactMode
                            figure(f4);
                            hold on;
                            plot(t90s(peakCounter), vsmooth(n90s(peakCounter)), 'r.', 'markersize', 10);
                            plot(t50s(peakCounter), vsmooth(n50s(peakCounter)), 'r.', 'markersize', 10);
                            plot(t10s(peakCounter), vsmooth(n10s(peakCounter)), 'r.', 'markersize', 10);
                            for i = 1:3
                                p5 = plot(peakTime, vsmooth(peakIndex), '^', 'markersize', 20, 'markerfacecolor', 'r', 'markeredgecolor', 'r');
                                pause(pauseTime);
                                delete(p5);
                                pause(pauseTime);
                            end
                            hold off;
                        end
                    end
                    % Coda: If the local maximum fails in any of the 5 tests, then it is rejected.
                    
                else
                    ampFull = vsmooth(peakIndex) - trough;
                    if peakCounter && ampFull > Amplobound && ampFull < Ampupbound && finalisedPeaks(peakCounter) > vsmooth(peakIndex) &&...
                            peakIndex - finalisedPeakIndices(peakCounter) <= 3*nSWupbound
                        
                        % 4/5. Test the rise time spread and check whether the baseline does not deviate too far from the 10% rise time mark:
                        if interactMode
                            [t1090, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, t(iSWstart - 1 + iTrough : peakIndex),...
                                vsmooth(iSWstart - 1 + iTrough : peakIndex), trough, interactMode, pauseTime);
                        else
                            [t1090, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, t(iSWstart - 1 + iTrough : peakIndex),...
                                vsmooth(iSWstart - 1 + iTrough : peakIndex), trough);
                        end
                        part1RT = i50 - i10;
                        part2RT = i90 - i50;
                        if part1RT/part2RT < 2.5 && part1RT/part2RT > .3 && i10 < nBLduration
                            peakCounter = peakCounter + 1;                  % A peak has been established
                            finalisedPeaks(peakCounter) = vsmooth(peakIndex);
                            peakTimes(peakCounter) = peakTime;
                            finalisedPeakIndices(peakCounter) = peakIndex;
                            amplitudes(peakCounter) = ampFull;
                            baselines(peakCounter) = trough;
                            tBaselines(peakCounter,:) = [t(iSWstart - 1 + iTrough) - dt, t(iSWstart - 1 + iTrough) + dt];
                            nBaselines(peakCounter,:) = [iSWstart - 1 + iTrough - 1, iSWstart - 1 + iTrough + 1];
                            elementRiseTimes(peakCounter) = peakIndex - iSWstart - 1 + iTrough;
                            riseTimes(peakCounter) = elementRiseTimes(peakCounter)*dt;
                            riseTimes1090(peakCounter) = t1090;
                            t10s(peakCounter) = t10;
                            t50s(peakCounter) = t50;
                            t90s(peakCounter) = t90;
                            n10s(peakCounter) = iSWstart - 1 + iTrough - 1 + i10;
                            n50s(peakCounter) = iSWstart - 1 + iTrough - 1 + i50;
                            n90s(peakCounter) = iSWstart - 1 + iTrough - 1 + i90;
                            v10s(peakCounter) = vsmooth(n10s(peakCounter));
                            v50s(peakCounter) = vsmooth(n50s(peakCounter));
                            v90s(peakCounter) = vsmooth(n90s(peakCounter));
                            decays(peakCounter) = calcDecay(vsmooth, baselines(peakCounter), peakIndex, dt);
                            
                            % PLOT3.1 - Bracket out the baseline:
                            if interactMode
                                figure(f4);
                                hold on;
                                delete([p3, p4]);
                                plot(tBaselines(peakCounter,1), trough, 'c.', 'markersize', 10);
                                plot(tBaselines(peakCounter,2), trough, 'b.', 'markersize', 10);
                                pause(.5*pauseTime);
                                
                                % PLOT4 - Indicate that the peak has been just established:
                                figure(f4);
                                plot(t90s(peakCounter), vsmooth(n90s(peakCounter)), 'r.', 'markersize', 10);
                                plot(t50s(peakCounter), vsmooth(n50s(peakCounter)), 'r.', 'markersize', 10);
                                plot(t10s(peakCounter), vsmooth(n10s(peakCounter)), 'r.', 'markersize', 10);
                                for i = 1:3
                                    p5 = plot(peakTime, vsmooth(peakIndex), '^', 'markersize', 20, 'markerfacecolor', 'r', 'markeredgecolor', 'r');
                                    pause(pauseTime);
                                    delete(p5);
                                    pause(pauseTime);
                                end
                                hold off;
                            end
                        end
                        % Coda: If the local maximum fails in any of the 5 tests, then it is rejected.
                        
                    elseif interactMode
                        % PLOT5 - Destroy unsuccessful mini candidate:
                        figure(f4);
                        hold on;
                        delete([p3, p4]);
                        p5 = plot(peakTime, vsmooth(peakIndex), 'x', 'markersize', 20, 'color', 'k', 'LineWidth', 2);
                        pause(2*pauseTime);
                        delete(p2);
                        delete(p5);
                        hold off;
                    end
                end
                
            else
                % PLOT5 - Destroy unsuccessful mini candidate:
                if interactMode
                    figure(f4);
                    hold on;
                    delete([p3, p4]);
                    p5 = plot(peakTime, vsmooth(peakIndex), 'x', 'markersize', 20, 'color', 'k', 'LineWidth', 2);
                    pause(2*pauseTime);
                    delete(p2);
                    delete(p5);
                    hold off;
                end
            end
        end
    end
end

% Assigning output data:
%            1               2          3                     4           5          6-7         8-9         10                11
coreMinis = [finalisedPeaks, peakTimes, finalisedPeakIndices, amplitudes, baselines, tBaselines, nBaselines, elementRiseTimes, riseTimes,...
    riseTimes1090, t10s, t50s, t90s, n10s, n50s, n90s, v10s, v50s, v90s, decays];
%   12             13    14    15    16    17    18    19    20    21    22
end





function coreMinis = detectionAlgorithmParallel(t, vsmooth, searchParameters, peakIndices, originalPeakInd)
% DETECTIONALGORITHMPARALLEL is a parallel minis detection algorithm and is
% a helper subfunction of DETECTMINIS.
%



% Assign major variables:
RTtype = searchParameters.RTinterval;                                       % The rise time type (20-80% or 10-90%)
dt = searchParameters.sampleInterval;                                       % The sample interval, ms
SWupbound = searchParameters.SWstart;                                       % The beginning of the search window, ms
nSWupbound = floor(SWupbound/dt);                                           % The element-wise beginning of the search window
tBLduration = searchParameters.BLduration;                                  % The baseline duration, ms
nBLduration = round(tBLduration/dt);                                        % The element-wise baseline duration
tRefractory = searchParameters.refractoryPeriod;                            % The refractory period, ms
nRefractory = ceil(tRefractory/dt);                                         % The element-wise refractory period
Amplobound = searchParameters.Amplobound;                                   % The lower bound on the amplitude size, mV
Ampupbound = searchParameters.Ampupbound;                                   % The upper bound on the amplitude size, mV
lengthTempStore = length(peakIndices);


% Initialise storage variables:
peakCounter = 0;                                                            % Established peak counter variable
finalisedPeaks = zeros(lengthTempStore,1);                                  % The array of confirmed peaks
peakTimes = zeros(lengthTempStore,1);
finalisedPeakIndices = zeros(lengthTempStore,1);                            % The indices of confirmed peaks
amplitudes = zeros(lengthTempStore,1);                                      % and their associated characteristics (below)
baselines = zeros(lengthTempStore,1);
tBaselines = zeros(lengthTempStore,2);
nBaselines = zeros(lengthTempStore,2);
elementRiseTimes = zeros(lengthTempStore,1);
riseTimes = zeros(lengthTempStore,1);
riseTimes1090 = zeros(lengthTempStore,1);
t10s = zeros(lengthTempStore,1);
t50s = zeros(lengthTempStore,1);
t90s = zeros(lengthTempStore,1);
n10s = zeros(lengthTempStore,1);
n50s = zeros(lengthTempStore,1);
n90s = zeros(lengthTempStore,1);
v10s = zeros(lengthTempStore,1);
v50s = zeros(lengthTempStore,1);
v90s = zeros(lengthTempStore,1);
decays = NaN(lengthTempStore,1);


% Major loop:
for iPeak = 1:length(peakIndices)
    peakIndex = peakIndices(iPeak);
    peakTime = t(peakIndex);
    
    % Filter out peaks that do not meat the inclusion criteria:
    % 1. Test whether the local maximum (peak) is the highest peak in the vicinity ahead:
    iEndRefractory = peakIndex + nRefractory;
    refractoryPeakIndices = peakIndices(peakIndices > peakIndex & peakIndices <= iEndRefractory);
    refractoryPeaks = vsmooth(refractoryPeakIndices);
    maxRefractoryPeak = max(refractoryPeaks);
    if isempty(maxRefractoryPeak) || vsmooth(peakIndex) > maxRefractoryPeak
        
        % 2. Test whether the local maximum (peak) has no established peaks in the vicinity behind:
        if peakCounter && peakTime - peakTimes(peakCounter) > tRefractory || ~peakCounter
            
            % 3. Test whether the local maximum (peak) is in the interval of acceptable amplitudes:
            % 3.1. Find the baseline of the peak:
            if peakCounter
                iSWstart = max([peakIndex - nSWupbound + 1, finalisedPeakIndices(peakCounter)]);
            else
                iSWstart = max([peakIndex - nSWupbound + 1, round(.5*nSWupbound)]);
            end
            [trough, iTrough] = min(vsmooth(iSWstart:peakIndex));
            iTrialBLend = min([iSWstart - 1 + iTrough + round(.2*nBLduration) peakIndex-4]);
            if peakCounter
                iTrialBLstart = max([iTrialBLend - nBLduration + 1 finalisedPeakIndices(peakCounter)]);
            else
                iTrialBLstart = iTrialBLend - nBLduration + 1;
            end
            trialBL = mean(vsmooth(iTrialBLstart:iTrialBLend));
            
            % 3.2. Test whether there are no higher peaks between the baseline and the proper peak (for removing edge effects):
            if isempty(vsmooth(originalPeakInd >= iSWstart - 1 + iTrough & originalPeakInd < peakIndex - 1))...
                    || ~(max(vsmooth(originalPeakInd(originalPeakInd >= iSWstart - 1 + iTrough & originalPeakInd < peakIndex - 1))) >= vsmooth(peakIndex))
                
                % 3.3. Test the amplitude:
                amp = vsmooth(peakIndex) - trialBL;
                if amp > Amplobound && amp < Ampupbound
                    
                    % 4. Test the rise time spread:
                    [t1090, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, t(iTrialBLend:peakIndex), vsmooth(iTrialBLend:peakIndex), trialBL);
                    part1RT = i50 - i10;
                    part2RT = i90 - i50;
                    
                    if part1RT/part2RT >= 5                                 % In case the first half is considerably larger than the second one
                        % 4.1. Correct the baseline if possible and re-test the amplitude:
                        newSWstart = peakIndices(find(peakIndices >= iTrialBLend, 1));
                        if newSWstart < peakIndex
                            [pass, amp, trialBL, iTrialBLstart, iTrialBLend, t1090, t10, i10, t50, i50, t90, i90] = correctBL(vsmooth, t, RTtype,...
                                newSWstart, peakIndex, nBLduration, Amplobound, Ampupbound, peakCounter, finalisedPeakIndices, false);
                            if pass
                                part1RT = i50 - i10;
                                part2RT = i90 - i50;
                                
                                if part1RT/part2RT >= 5                     % In case the first half is still considerably larger than the second one
                                    % 4.2. Repeat the correction:
                                    newSWstart = peakIndices(find(peakIndices >= iTrialBLend, 1));
                                    if newSWstart ~= peakIndex
                                        [pass, amp, trialBL, iTrialBLstart, iTrialBLend, t1090, t10, i10, t50, i50, t90, i90] = correctBL(vsmooth, t,...
                                            RTtype, newSWstart, peakIndex, nBLduration, Amplobound, Ampupbound, peakCounter, finalisedPeakIndices, false);
                                    end
                                    if ~pass
                                        continue
                                    end
                                end    
                            else
                                continue
                            end
                        end
                        
                    elseif part1RT/part2RT <= .2                            % In case the first half is considerably smaller than the second one
                        % 4.3. Shift the peak backwards if possible:
                        newPeakIndex = peakIndices(find(peakIndices < peakIndex, 1, 'last'));
                        if ~isempty(newPeakIndex) && newPeakIndex > iTrialBLend
                            newPeakTime = t(newPeakIndex);
                            
                            % 4.4. Re-test the amplitude:
                            newAmp = vsmooth(newPeakIndex) - trialBL;
                            if newAmp > Amplobound && newAmp < Ampupbound
                                peakIndex = newPeakIndex;
                                peakTime = newPeakTime;
                                amp = newAmp;
                                [t1090, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, t(iTrialBLend:peakIndex), vsmooth(iTrialBLend:peakIndex), trialBL);
                            end
                        end
                    end
                    
                    % 5. Check whether the baseline does not deviate too far from the 10% rise time mark:
                    if (strcmpi(RTtype, '10-90%') && iTrialBLstart - 1 + i10 - iTrialBLend > ceil(.5*nBLduration))...
                            || (strcmpi(RTtype, '20-80%') && iTrialBLstart - 1 + i10 - iTrialBLend > nBLduration)
                        
                        % 5.1. Correct the baseline:
                        iSWstart = iTrialBLstart - 1 + i10 - ceil(.5*nBLduration);
                        [~, iTrough] = min(vsmooth(iSWstart:peakIndex));
                        iTrialBLend = min([iSWstart - 1 + iTrough + round(.2*nBLduration) peakIndex-4]);
                        iTrialBLstart = iTrialBLend - nBLduration + 1;
                        trialBL = mean(vsmooth(iTrialBLstart:iTrialBLend));
                        
                        % 5.2. Test the new amplitude:
                        amp = vsmooth(peakIndex) - trialBL;
                        if amp > Amplobound && amp < Ampupbound
                            [t1090, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, t(iTrialBLend:peakIndex), vsmooth(iTrialBLend:peakIndex), trialBL);
                            peakCounter = peakCounter + 1;                  % A peak has been established
                            finalisedPeaks(peakCounter) = vsmooth(peakIndex);
                            peakTimes(peakCounter) = peakTime;
                            finalisedPeakIndices(peakCounter) = peakIndex;
                            amplitudes(peakCounter) = amp;
                            baselines(peakCounter) = trialBL;
                            tBaselines(peakCounter,:) = [t(iTrialBLstart), t(iTrialBLend)];
                            nBaselines(peakCounter,:) = [iTrialBLstart, iTrialBLend];
                            elementRiseTimes(peakCounter) = peakIndex - iTrialBLstart + 1;
                            riseTimes(peakCounter) = elementRiseTimes(peakCounter)*dt;
                            riseTimes1090(peakCounter) = t1090;
                            t10s(peakCounter) = t10;
                            t50s(peakCounter) = t50;
                            t90s(peakCounter) = t90;
                            n10s(peakCounter) = iTrialBLend - 1 + i10;
                            n50s(peakCounter) = iTrialBLend - 1 + i50;
                            n90s(peakCounter) = iTrialBLend - 1 + i90;
                            v10s(peakCounter) = vsmooth(n10s(peakCounter));
                            v50s(peakCounter) = vsmooth(n50s(peakCounter));
                            v90s(peakCounter) = vsmooth(n90s(peakCounter));
                            decays(peakCounter) = calcDecay(vsmooth, baselines(peakCounter), peakIndex, dt);
                        end
                        % Coda: If the local maximum fails in any of the 5 tests, then it is rejected.
                        
                    else
                        peakCounter = peakCounter + 1;                      % A peak has been established
                        finalisedPeaks(peakCounter) = vsmooth(peakIndex);
                        peakTimes(peakCounter) = peakTime;
                        finalisedPeakIndices(peakCounter) = peakIndex;
                        amplitudes(peakCounter) = amp;
                        baselines(peakCounter) = trialBL;
                        tBaselines(peakCounter,:) = [t(iTrialBLstart), t(iTrialBLend)];
                        nBaselines(peakCounter,:) = [iTrialBLstart, iTrialBLend];
                        elementRiseTimes(peakCounter) = peakIndex - iTrialBLstart + 1;
                        riseTimes(peakCounter) = elementRiseTimes(peakCounter)*dt;
                        riseTimes1090(peakCounter) = t1090;
                        t10s(peakCounter) = t10;
                        t50s(peakCounter) = t50;
                        t90s(peakCounter) = t90;
                        n10s(peakCounter) = iTrialBLend - 1 + i10;
                        n50s(peakCounter) = iTrialBLend - 1 + i50;
                        n90s(peakCounter) = iTrialBLend - 1 + i90;
                        v10s(peakCounter) = vsmooth(n10s(peakCounter));
                        v50s(peakCounter) = vsmooth(n50s(peakCounter));
                        v90s(peakCounter) = vsmooth(n90s(peakCounter));
                        decays(peakCounter) = calcDecay(vsmooth, baselines(peakCounter), peakIndex, dt);
                    end
                    % Coda: If the local maximum fails in any of the 5 tests, then it is rejected.
                    
                else
                    ampFull = vsmooth(peakIndex) - trough;
                    if peakCounter && ampFull > Amplobound && ampFull < Ampupbound && finalisedPeaks(peakCounter) > vsmooth(peakIndex) &&...
                            peakIndex - finalisedPeakIndices(peakCounter) <= 3*nSWupbound
                        
                        % 4/5. Test the rise time spread and check whether the baseline does not deviate too far from the 10% rise time mark:
                        [t1090, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, t(iSWstart - 1 + iTrough : peakIndex),...
                            vsmooth(iSWstart - 1 + iTrough : peakIndex), trough);
                        part1RT = i50 - i10;
                        part2RT = i90 - i50;
                        if part1RT/part2RT < 2.5 && part1RT/part2RT > .3 && i10 < nBLduration
                            peakCounter = peakCounter + 1;                  % A peak has been established
                            finalisedPeaks(peakCounter) = vsmooth(peakIndex);
                            peakTimes(peakCounter) = peakTime;
                            finalisedPeakIndices(peakCounter) = peakIndex;
                            amplitudes(peakCounter) = ampFull;
                            baselines(peakCounter) = trough;
                            tBaselines(peakCounter,:) = [t(iSWstart - 1 + iTrough) - dt, t(iSWstart - 1 + iTrough) + dt];
                            nBaselines(peakCounter,:) = [iSWstart - 1 + iTrough - 1, iSWstart - 1 + iTrough + 1];
                            elementRiseTimes(peakCounter) = peakIndex - iSWstart - 1 + iTrough;
                            riseTimes(peakCounter) = elementRiseTimes(peakCounter)*dt;
                            riseTimes1090(peakCounter) = t1090;
                            t10s(peakCounter) = t10;
                            t50s(peakCounter) = t50;
                            t90s(peakCounter) = t90;
                            n10s(peakCounter) = iSWstart - 1 + iTrough - 1 + i10;
                            n50s(peakCounter) = iSWstart - 1 + iTrough - 1 + i50;
                            n90s(peakCounter) = iSWstart - 1 + iTrough - 1 + i90;
                            v10s(peakCounter) = vsmooth(n10s(peakCounter));
                            v50s(peakCounter) = vsmooth(n50s(peakCounter));
                            v90s(peakCounter) = vsmooth(n90s(peakCounter));
                            decays(peakCounter) = calcDecay(vsmooth, baselines(peakCounter), peakIndex, dt);
                        end
                        % Coda: If the local maximum fails in any of the 5 tests, then it is rejected.
                        
                    end
                end
            end
        end
    end
end

% Assigning output data:
%            1               2          3                     4           5          6-7         8-9         10                11
coreMinis = [finalisedPeaks, peakTimes, finalisedPeakIndices, amplitudes, baselines, tBaselines, nBaselines, elementRiseTimes, riseTimes,...
    riseTimes1090, t10s, t50s, t90s, n10s, n50s, n90s, v10s, v50s, v90s, decays];
%   12             13    14    15    16    17    18    19    20    21    22
end





function [minis, f, p] = manualAdjust(minis, t, V, f, p, figureAxes, tBLduration, RTtype)
% MANUALADJUST is a helper subfunction of DETECTMINIS. It allows the user
% to manually edit the events detected by the detectMinis program.
%
%   [MINIS, F, P] = MANUALADJUST(minis, t, V, f, p, figureAxes, tBLduration)
%   controls the manual event edition. MINIS is a matrix output of the
%   minis detection algorithm as used by the detectMinis program. T is a
%   time vector (ms). V is electrophysiological recording vector (mV or
%   nA). F is the detected events summary figure handle. P contains handles
%   to the graphical objects in the figure. FIGUREAXES is the axes handle
%   of the summary figure. TBLDURATION is a scalar variable representing
%   the baseline duration in miliseconds (ms). RTTYPE is the rise-time
%   type: '10-90%' or '20-80%'.
%   The output variables MINIS correspond to the input variable but after
%   the manual event edition. F is the summary figure handle with P being
%   the set of handles of the graphical objects in the figure.
%
%   The events are edited using left and right mouse buttons and certain
%   keyboard buttons:
%   'Left mouse key' - Add an event or correct the baseline of an event.
%       Left-click near the peak you want to be detected as a minis-like
%       event. In this case, the program will automatically find the
%       baseline position. Otherwise, if you want to edit the baseline,
%       left-click between the baseline marks. The Baseline marks would
%       turn into magenta colour. Then click on your chosen through to
%       position the new baseline.
%   'Right mouse key' - Delete an event. Right-click near the peak of the
%       mini that you want to be deleted.
%   'Carriage return (enter)' and 'Esc' exits the manual editing mode. All
%       changes are saved and returned.
%   'Right arrow key' - Next interval.
%   'Left arrow key' - Previous interval.
%   'Up arrow key' - zoom out (also 'o' key).
%   'Down arrow key' - zoom in (also 'i' key).
%   's' - Return to the very beginning (the first interval).
%   'e' - Go to the very end (the last interval).
%   'i' - Zoom in.
%   'o' - Zoom out.
%   'h' - Bring back the help window.
%



% Set axis limits:
XLim = [t(1) t(20000)];
figure(f);
hold on
set(figureAxes, 'XLim', XLim);
zoomFactor = .3;
p2 = zeros(1,2);

dt = t(2)-t(1);
nBLduration = round(tBLduration/dt);

% Edit minis:
edit = 1;
while edit
    try
        [x, y, button] = ginput(1);
    catch %#ok<CTCH>
    end
    if isempty(button) || button == 27
        XLim = [t(1) t(end)];
        set(figureAxes, 'XLim', XLim);
        return
    end
    XLim = get(figureAxes, 'XLim');
    YLim = get(figureAxes, 'YLim');
    if (button == 1 || button == 3) && (x < XLim(1) || x > XLim(2) || y < YLim(1) || y > YLim(2))
        button = 0;
    end
    ix = round(x/dt);
    switch button
        case 1 % Add a mini or correct the baseline
            % Recognise whether it is a baseline correction or a mini:
            iUpBL = find(minis(:,9) >= ix, 1);
            if minis(iUpBL,8) <= ix
                % Mark the baseline for correction:
                figure(f);
                p2(1) = plot(minis(iUpBL,6), minis(iUpBL,5), 'm.', 'markersize', 10);
                p2(2) = plot(minis(iUpBL,7), minis(iUpBL,5), 'm.', 'markersize', 10);
                % Correct the baseline:
                editBL = 1;
                while editBL
                    [x, ~, button] = ginput(1);
                    if button == 1
                        ix = round(x/dt);
                        editBL = 0;
                    end
                end
                iBLend = max([min([ix - 1 + round(.5*nBLduration) length(t)]) nBLduration]);
                iBLstart = iBLend - nBLduration + 1;
                BL = mean(V(iBLstart:iBLend));
                if minis(iUpBL,1) - BL > 0
                    % Estimate the new rise times:
                    [t1090, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, t(iBLstart:minis(iUpBL,3)), V(iBLstart:minis(iUpBL,3)), BL);
                    minis(iUpBL,:) = [minis(iUpBL,1:3), minis(iUpBL,1)-BL, BL, t(iBLstart), t(iBLend), iBLstart, iBLend, minis(iUpBL,3)-iBLstart+1,...
                        (minis(iUpBL,3)-iBLstart+1)*dt, t1090, t10, t50, t90, iBLstart-1+i10, iBLstart-1+i50, iBLstart-1+i90, V(iBLstart-1+i10),...
                        V(iBLstart-1+i50), V(iBLstart-1+i90), minis(iUpBL,22)];
                end
                delete(p2);
            else
                % Localise a mini:
                iUp = find(minis(:,3) >= ix, 1);
                iPeakUpbound = minis(iUp,3);
                if isempty(iPeakUpbound)
                    iPeakUpbound = length(t);
                end
                iLo = find(minis(:,3) <= ix, 1, 'last');
                iPeakLobound = minis(iLo,3);
                prevPeak = 1;
                if isempty(iPeakLobound)
                    iPeakLobound = 1;
                    prevPeak = 0;
                end
                [peaks, locs] = findpeaks(double(V(iPeakLobound:iPeakUpbound)));
                if isempty(peaks)
                    continue % do nothing
                end
                locs = iPeakLobound + locs - 1;
                distances = abs(locs - ix);
                [~, iNearestPeak] = min(distances);
                peak = peaks(iNearestPeak);
                iPeak = locs(iNearestPeak);
                if iPeak ~= iPeakLobound && iPeak ~= iPeakUpbound
                    % Measure the baseline:
                    [~, iTrough] = min(V(iPeakLobound:iPeak));
                    iBLend = min([iPeakLobound - 1 + iTrough iPeak]);
                    if prevPeak
                        iBLstart = max([iBLend - nBLduration + 1 iPeakLobound]);
                    else
                        iBLstart = max([iBLend - nBLduration + 1 1]);
                    end
                    BL = mean(V(iBLstart:iBLend));
                    if peak - BL <= 0
                        continue
                    end
                    % Estimate the rise times:
                    [t1090, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, t(iBLend:iPeak), V(iBLend:iPeak), BL);
                    % Test whether the baseline is at the foothill:
                    if iBLstart - 1 + i10 - iBLend > ceil(.75*nBLduration)
                        iPeakLobound = iBLstart - 1 + i10 - ceil(.75*nBLduration);
                        [~, iTrough] = min(V(iPeakLobound:iPeak));
                        iBLend = min([iPeakLobound - 1 + iTrough iPeak]);
                        iBLstart = iBLend - nBLduration + 1;
                        BL = mean(V(iBLstart:iBLend));
                        [t1090, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, t(iBLend:iPeak), V(iBLend:iPeak), BL);
                    end
                    decay = calcDecay(V, BL, iPeak, dt);
                    % Establish the new mini:
                    %          1               2          3                     4           5          6-7         8-9         10                11
                    %         [finalisedPeaks, peakTimes, finalisedPeakIndices, amplitudes, baselines, tBaselines, nBaselines, elementRiseTimes, riseTimes,...
                    %   riseTimes1090, t10s, t50s, t90s, n10s, n50s, n90s, v10s, v50s, v90s, decay];
                    %   12             13    14    15    16    17    18    19    20    21    22
                    
                    %          1     2         3      4        5   6            7          8         9       10                11
                    newData = [peak, t(iPeak), iPeak, peak-BL, BL, t(iBLstart), t(iBLend), iBLstart, iBLend, iPeak-iBLstart+1, (iPeak-iBLstart+1)*dt,...
                        t1090, t10, t50, t90, iBLend-1+i10, iBLend-1+i50, iBLend-1+i90, V(iBLend-1+i10), V(iBLend-1+i50), V(iBLend-1+i90), decay];
                    %   12     13   14   15   16            17            18            19               20               21               22
                    if isempty(iLo)
                        minis = [newData; minis]; %#ok<AGROW>
                    elseif isempty(iUp)
                        minis = [minis; newData]; %#ok<AGROW>
                    else
                        minis = [minis(1:iLo,:); newData; minis(iUp:end,:)];
                    end
                end
            end
            % Re-plot minis:
            delete(p);
            p(1) = plot(minis(:,2), minis(:,1), '^', 'markersize', 10, 'markeredgecolor', 'r');
            p(2) = plot(minis(:,15), minis(:,21), 'r.', 'markersize', 10);
            p(3) = plot(minis(:,14), minis(:,20), 'r.', 'markersize', 10);
            p(4) = plot(minis(:,13), minis(:,19), 'r.', 'markersize', 10);
            p(5) = plot(minis(:,6), minis(:,5), 'c.', 'markersize', 10);
            p(6) = plot(minis(:,7), minis(:,5), 'b.', 'markersize', 10);
            if exist('fixedYLim', 'var')
              set(figureAxes, 'YLim', fixedYLim);
            end
            
        case 3 % Delete a mini
            distances = abs(minis(:,3) - ix);
            [~, miniToDelete] = min(distances);
            minis(miniToDelete,:) = [];
            delete(p);
            p(1) = plot(minis(:,2), minis(:,1), '^', 'markersize', 10, 'markeredgecolor', 'r');
            p(2) = plot(minis(:,15), minis(:,21), 'r.', 'markersize', 10);
            p(3) = plot(minis(:,14), minis(:,20), 'r.', 'markersize', 10);
            p(4) = plot(minis(:,13), minis(:,19), 'r.', 'markersize', 10);
            p(5) = plot(minis(:,6), minis(:,5), 'c.', 'markersize', 10);
            p(6) = plot(minis(:,7), minis(:,5), 'b.', 'markersize', 10);
            if exist('fixedYLim', 'var')
              set(figureAxes, 'YLim', fixedYLim);
            end
            
        case 29 % Next window (right arrow key)
            prevXLim = get(figureAxes,'XLim');
            interval = prevXLim(2) - prevXLim(1);
            XLim = [prevXLim(2)-.1*interval prevXLim(2)+.9*interval];
            if XLim(2) > t(end)
                XLim = [t(end)-interval t(end)];
            end
            set(figureAxes, 'XLim', XLim);
            if exist('fixedYLim', 'var')
              set(figureAxes, 'YLim', fixedYLim);
            end
            
        case 28 % Previous window (left arrow key)
            currXLim = get(figureAxes,'XLim');
            interval = currXLim(2) - currXLim(1);
            XLim = [currXLim(1)-.9*interval currXLim(1)+.1*interval];
            if XLim(1) < t(1)
                XLim = [t(1) interval];
            end
            set(figureAxes, 'XLim', XLim);
            if exist('fixedYLim', 'var')
              set(figureAxes, 'YLim', fixedYLim);
            end
            
        case 's' % Go to the beginning
            currXLim = get(figureAxes,'XLim');
            XLim = [t(1) currXLim(2)-currXLim(1)];
            set(figureAxes, 'XLim', XLim);
            if exist('fixedYLim', 'var')
              set(figureAxes, 'YLim', fixedYLim);
            end
            
        case 'e' % Go to the end
            prevXLim = get(figureAxes,'XLim');
            XLim = [t(end)-prevXLim(2)-prevXLim(1) t(end)];
            set(figureAxes, 'XLim', XLim);
            if exist('fixedYLim', 'var')
              set(figureAxes, 'YLim', fixedYLim);
            end

        case 'f' % Fix y-axis limits
            if exist('fixedYLim', 'var')
              clear fixedYLim
              ylim("auto");
            else
              fixedYLim = get(figureAxes,'YLim');
            end
            
        case {'i', 31} % Zoom in
            prevXLim = get(figureAxes,'XLim');
            XLim = [prevXLim(1)+zoomFactor*(prevXLim(2)-prevXLim(1)) prevXLim(2)-zoomFactor*(prevXLim(2)-prevXLim(1))];
            set(figureAxes, 'XLim', XLim);
            if exist('fixedYLim', 'var')
              set(figureAxes, 'YLim', fixedYLim);
            end
            
        case {'o', 30} % Zoom out
            prevXLim = get(figureAxes,'XLim');
            XLim = [prevXLim(1)-zoomFactor*(prevXLim(2)-prevXLim(1)) prevXLim(2)+zoomFactor*(prevXLim(2)-prevXLim(1))];
            set(figureAxes, 'XLim', XLim);
            if exist('fixedYLim', 'var')
              set(figureAxes, 'YLim', fixedYLim);
            end
            
        case 'h' % Help
            bringHelp;
    end
end
end





function bringHelp
% BRINGHELP is a helper subfunction of DETECTMINIS. It contains the help
% text for the manual adjustment of detected minis-like events using the
% mouse and the keyboard.
%



helpText1 =  'The events are edited using left and right mouse buttons and certain';
helpText2 =  'keyboard buttons:';
helpText3 =  ' ';
helpText4 =  'Left mouse key  -  Add an event. Left-click near the peak you want to';
helpText5 =  '     be detected as a mini. The program will automatically find the';
helpText6 =  '     baseline position. If you want to edit the baseline, left-click';
helpText7 =  '     between the baseline marks. The Baseline marks would turn into';
helpText8 =  '     magenta colour. Then click on your chosen through to position the';
helpText9 =  '     new baseline.';
helpText10 = 'Right mouse key  -  Delete an event. Right-click near the peak of the';
helpText11 = '     mini that you want to be deleted.';
helpText12 = 'Carriage return (enter) and Esc exits the manual editing mode. All';
helpText13 = '     changes are saved and returned.';
helpText14 = 'Right arrow key  -  Next interval.';
helpText15 = 'Left arrow key  -  Previous interval.';
helpText16 = 'Up arrow key - Zoom out (also o key).';
helpText17 = 'Down arrow key - Zoom in (also i key).';
helpText18 = 's  -  Return to the very beginning (the first interval).';
helpText19 = 'e  -  Go to the very end (the last interval).';
helpText20 = 'i or downward arrow key  -  Zoom in.';
helpText21 = 'o or upward arrow key -  Zoom out.';
helpText22 = 'f  -  Freeze/unfreeze y-axis limits.';
helpText23 = 'h  -  Bring back the help window.';
helpText = char(helpText1,helpText2,helpText3,helpText4,helpText5,helpText6,helpText7,helpText8,helpText9,helpText10,helpText11,...
    helpText12,helpText13,helpText14,helpText15,helpText16,helpText17,helpText18,helpText19,helpText20,helpText21,helpText22,helpText23);
uiwait(helpdlg(helpText,'Help for Manual Event Editing'));
end





function [pass, amp, trialBL, iTrialBLstart, iTrialBLend, t1090, t10, i10, t50, i50, t90, i90, varargout] = correctBL(vsmooth, t, RTtype, newSWstart,...
    peakIndex, nBLduration, Amplobound, Ampupbound, peakCounter, finalisedPeakIndices, interactMode, varargin)
% CORRECTBL is a helper subfunction of DETECTMINIS.
%



iSWstart = newSWstart+1;
[~, iTrough] = min(vsmooth(iSWstart:peakIndex));
iTrialBLend = min([iSWstart - 1 + iTrough + round(.2*nBLduration) peakIndex-4]);
if peakCounter
    iTrialBLstart = max([iTrialBLend - nBLduration + 1 finalisedPeakIndices(peakCounter)]);
else
    iTrialBLstart = iTrialBLend - nBLduration + 1;
end
trialBL = mean(vsmooth(iTrialBLstart:iTrialBLend));

% PLOT3.2 - Bracket out the new baseline:
if interactMode
    pauseTime = cell2mat(varargin(1));
    f4 = cell2mat(varargin(2));
    p3 = cell2mat(varargin(3));
    p4 = cell2mat(varargin(4));
    figure(f4);
    hold on;
    delete([p3 p4]);
    p3 = plot(t(iTrialBLstart), trialBL, 'c.', 'markersize', 10);
    p4 = plot(t(iTrialBLend), trialBL, 'b.', 'markersize', 10);
    hold off;
    pause(.5*pauseTime);
    varargout{1} = {p3};
    varargout{2} = {p4};
end

% 6.2. Re-test the amplitude:
amp = vsmooth(peakIndex) - trialBL;
if amp > Amplobound && amp < Ampupbound
    pass = true;
    if interactMode
        [t1090, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, t(iTrialBLend:peakIndex), vsmooth(iTrialBLend:peakIndex), trialBL, interactMode, pauseTime);
    else
        [t1090, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, t(iTrialBLend:peakIndex), vsmooth(iTrialBLend:peakIndex), trialBL);
    end
else
    pass = false;
    t1090 = [];
    t10 = [];
    i10 = [];
    t50 = [];
    i50 = [];
    t90 = [];
    i90 = [];
end
end





function decay = calcDecay(V, baseline, peakIndex, dt)
% CALCDECAY is a helper subfunction of DETECTMINIS.
%



V = V - baseline;
maxDecay = 30; % milliseconds
maxDecay = maxDecay/dt;
endMaxDecay = min([peakIndex+maxDecay numel(V)]);
vsmoothOI = V(peakIndex:endMaxDecay);
iDecay = find(vsmoothOI <= vsmoothOI(1)/exp(1), 1);
if isempty(iDecay)
    decay = NaN;
else
    decay = iDecay*dt;
end
end