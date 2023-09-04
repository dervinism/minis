function manualFittingMinis(wd, targetFilename, noiseFilename, targetExcludedTimes, noiseExcludedTimes, detectionParameters, simulationParameters, classificationParameters, graphicsFormats, filtering, parallelCores)
% MANUALFITTINGMINIS provides a graphical user interface for manually
% finding the distribution of miniature excitatory postsynaptic potentials
% (mEPSPs) or minis (or minis-like events) given the separate recording
% traces of minis + noise and noise membrane potential (or current clamp)
% data.
%
%   MANUALFITTINGMINIS(wd, targetFilename, noiseFilename,
%   targetExcludedTimes, noiseExcludedTimes, detectionParameters,
%   simulationParameters, classificationParameters, graphicsFormats,
%   parallelCores) provides a GUI for the manual control of minis
%   distribution parameters. The minis are simulated using this
%   distribution and written onto the noise data in order to match it with
%   the noise + minis recording data. WD is a string variable containing
%   the working directory path. TARGETFILENAME is a string variable
%   containing a full name of the target (minis + noise) data file.
%   NOISEFILENAME is a string variable containing a full name of the noise
%   data file. TARGETEXCLUDEDTIMES is a vector with time values that should
%   be excluded from the minis detection analysis of the target data.
%   NOISEEXCLUDEDTIMES is a vector with time values that should be excluded
%   from the minis detection analysis of the noise data.
%   DETECTIONPARAMETERS is a structure variable with parameters used in the
%   minis detection analysis. For the description of the fields of
%   detectionParameters please see the detectMinis function.
%   SIMULATIONPARAMETERS is a structure variable containing parameters used
%   for simulating minis. For the description of simulationParameters
%   fields please see the simulateMinis function (corresponds to PARAMETERS
%   input variable). CLASSIFICATIONPARAMETERS is a structure variable
%   containing classification ranges for amplitudes and rise times. The
%   fields of this variable are
%       'amplitudeArray' - a vector of classification amplitudes;
%       'amplitudeArrayExt' - a vector of classification amplitudes with an
%           additional bin at the end;
%       'riseTimeArray' - a vector of classification rise times;
%       'riseTimeArrayExt' - a vector of classification rise times with an
%           additional bin at the end;
%   GRAPHICSFORMATS is a cell array of file formats with descriptions for
%   saving figures. FILTERING is a string variable that controls the
%   band-stop filtering of the electrophysiological data with values 'on'
%   and 'off' corresponding filtering and leaving the data unfiltered,
%   respectively.PARALLELCORES is a scalar variable with values set to
%   correspond to the number of processor cores intended to use for the
%   parallel processing in the minis detection task.
%

F = zeros(1,4);
waveformT.estimate = 1;
waveformT.riseTimeArray = classificationParameters.riseTimeArray;
waveformT.tau_m = simulationParameters.tau_m;
waveformT.classificationParameters = classificationParameters;
waveformN.estimate = false;
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

%% Detect events:
[minisTarget, ~, targetWave, targetProperties, F0] = detectMinisStandalone('off', targetFilename, detectionParameters,...
    targetExcludedTimes, classificationParameters, parallelCores, 'none', 0, waveformT, targetFilt);
tau_m = ones(size(minisTarget,1),1)*targetWave.parameters.tau_m;
minisTarget = [minisTarget tau_m];
[minisNoise, filtNoise, ~, noiseProperties, F(1:2)] = detectMinisStandalone('off', noiseFilename, detectionParameters,...
    noiseExcludedTimes, classificationParameters, parallelCores, 'none', 0, waveformN, noiseFilt);
noiseExcludedTimes = noiseProperties.excludedTimes;
noiseProperties.sweep = filtNoise;
clear filtering waveformN filtNoise


%% Save the event text files:
button = questdlg('Save the event text files?','Save File','Yes','No','Yes');
if strcmpi(button, 'Yes')
    % Target events:
    [filename, pathname, filterIndex] = uiputfile({'*.txt','Text files (*.txt)'},'Save Target Event Log as', 'detected target events.txt');
    if filterIndex
        fid = fopen(fullfile(pathname,filename),'wt+');
        fprintf(fid, '%16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\n',...
            'Peak potential', 'Peak Time', 'Peak Index', 'Amplitude', 'Baseline (BL)', 'BL start time', 'BL end time', 'BL start index', 'BL end index',...
            'Rise time (RT) length', 'RT', '10-90% RT', '10% RT time mark', '50% RT time mark', '90% RT time mark',...
            '10% RT index', '50% RT index', '90% RT index', '10% RT potential', '50% RT potential', '90% RT potential', 'tau');
        fprintf(fid,'%16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\n',minisTarget');
        fclose(fid);
    end
    % Noise events:
    [filename, pathname, filterIndex] = uiputfile({'*.txt','Text files (*.txt)'},'Save Noise Event Log as', 'detected noise events.txt');
    if filterIndex
        fid = fopen(fullfile(pathname,filename),'wt+');
        fprintf(fid, '%16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\t %16s\n',...
            'Peak potential', 'Peak Time', 'Peak Index', 'Amplitude', 'Baseline (BL)', 'BL start time', 'BL end time', 'BL start index', 'BL end index',...
            'Rise time (RT) length', 'RT', '10-90% RT', '10% RT time mark', '50% RT time mark', '90% RT time mark',...
            '10% RT index', '50% RT index', '90% RT index', '10% RT potential', '50% RT potential', '90% RT potential');
        fprintf(fid,'%16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\t %16.14g\n',minisNoise');
        fclose(fid);
    end
end


%% Save the detected event figure files:
button = questdlg('Save the graphs showing detected events?','Save File','Yes','No','Yes');
if strcmpi(button, 'Yes')
    [filename, pathname, filterIndex] = uiputfile(graphicsFormats,'Save Target Event Graph as', 'detected target events');
    if filterIndex
        figFullName = fullfile(pathname, filename);
        saveas(F0(1), figFullName);
    end
    [filename, pathname, filterIndex] = uiputfile(graphicsFormats,'Save Noise Event Graph as', 'detected noise events');
    if filterIndex
        figFullName = fullfile(pathname, filename);
        saveas(F(1), figFullName);
    end
end


%% Save the frequency spectrum graphs:
if numel(F0)>1 && numel(F)>1
    button = questdlg('Save the frequency spectrum graphs?','Save File','Yes','No','Yes');
    if strcmpi(button, 'Yes')
        [figureFilename, figurePathname, filterIndex] = uiputfile(graphicsFormats,'Save Target Frequency Spectrum as', 'target frequency spectrum');
        if filterIndex
            figFullName = fullfile(figurePathname, figureFilename);
            saveas(F0(2), figFullName);
        end
        [figureFilename, figurePathname, filterIndex] = uiputfile(graphicsFormats,'Save Noise Frequency Spectrum as', 'noise frequency spectrum');
        if filterIndex
            figFullName = fullfile(figurePathname, figureFilename);
            saveas(F(2), figFullName);
        end
    end
end


%% Close the graphs:
button = questdlg('Close the graphs?','Close window','Yes','No','Yes');
if strcmpi(button, 'Yes')
    if numel(F0)>1 && numel(F)>1
        close([F0 F(1:2)]);
    else
        close([F0(1) F(1)]);
    end
    F(1) = 0;
end

%% Initialise range variables:
amplitudeArray = classificationParameters.amplitudeArray;
amplitudeArrayExt = classificationParameters.amplitudeArrayExt;
riseTimeArray = classificationParameters.riseTimeArray;
riseTimeArrayExt = classificationParameters.riseTimeArrayExt;
lengthRatio = (length(targetProperties.sweep)-length(targetExcludedTimes))/(length(noiseProperties.sweep)-length(noiseExcludedTimes));


%% Initialise distribution variables and plot them:
[targetAmplitudes, ~, target2D] = classifyMinis(minisTarget(:,4), minisTarget(:,12), classificationParameters);
noiseAmplitudes = classifyMinis(minisNoise(:,4), minisNoise(:,12), classificationParameters);
minis1D = round(lengthRatio*(targetAmplitudes - noiseAmplitudes));
minis1D(minis1D < 0) = 0;
simulationParameters.tau_sy1 = 1;
simulationParameters.tau_sy2 = .1*simulationParameters.tau_sy1;
try
    simulationParameters = changeParametersGUI(simulationParameters);
catch err
    disp(err);
end
[~, V, ~, ~, shapes] = simulateMinis('Zero', simulationParameters, 1, 'Arbitrary', minis1D, noiseProperties, 'basal', parallelCores, amplitudeArray);
[~, ~, minis2D] = classifyMinis(shapes(:,2), shapes(:,3), classificationParameters);


%% Obtain first estimate based on subtracted distribution:
detectionParameters.sampleInterval = noiseProperties.dt;
detectionParameters.smoothWindow = round(detectionParameters.smoothWindow/detectionParameters.sampleInterval);
detectionParameters.smoothWindowLite = 8;
filtering.state = 'off';
%Ask whether to display the detection summary plot:
summaryBtn = questdlg('Do you want to display the detection summary plot (choosing NO saves memory)?','Figure display','Yes','No','No');
if strcmpi(summaryBtn,'Yes')
    options.summaryPlot = 'on';
    [minisSimulated, ~, ~, simWave, F(1)] = detectMinis(V, noiseExcludedTimes, detectionParameters, filtering, waveformT, parallelCores, options);
else
    options.summaryPlot = 'off';
    [minisSimulated, ~, ~, simWave] = detectMinis(V, noiseExcludedTimes, detectionParameters, filtering, waveformT, parallelCores, options);
end
tau_m = ones(size(minisSimulated,1),1)*simWave.parameters.tau_m;
minisSimulated = [minisSimulated tau_m];
simulated2D = lengthRatio*hist2d(minisSimulated(:,4), minisSimulated(:,12), amplitudeArrayExt, riseTimeArrayExt);
simulated2D = simulated2D(1:end-1,1:end-1);
[F(2), F(3), F(4)] = plotMinis(amplitudeArray, riseTimeArray, target2D, simulated2D, minis2D, shapes);


%% Open iteration GUI:
clear F0 amplitudeArray amplitudeArrayExt button minis1D minisNoise minisTarget noiseAmplitudes noiseExcludedTimes noiseFilename riseTimeArray...
    riseTimeArrayExt simWave summaryBtn targetAmplitudes targetExcludedTimes targetProperties targetWave tau_m waveformT
manualFitMinisGUI({wd}, V, minisSimulated, target2D, simulated2D, minis2D, lengthRatio, noiseProperties, detectionParameters,...
    simulationParameters, classificationParameters, graphicsFormats, parallelCores, F, options, shapes)