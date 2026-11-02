function [waveform, averageAmp, medianAmp, waveformRT, parameters, varargout] = averageWaveEffect(minis, V, dt, tau_m, excludedT,...
  tBLduration, RTtype, riseTimeArrayExt, draw, varargin)
% AVERAGEWAVEEFFECT is a helper function of DETECTMINIS. By default it
% estimates the average waveform using ALL detected minis-like events; pass
% TOP10PERCENTONLY (see below) to restrict it to the top 10% largest
% instead, as the preprocessing pipeline (dataMinis.m) does.
%
%   [WAVEFORM WAVEFORMAMP WAVEFORMRT PARAMETERS] = AVERAGEWAVEEFFECT(...
%       minis, V, dt, tau_m, excludedT, tBLduration, RTtype,...
%       riseTimeArray, draw)
%   Estimates the average waveform and the rise time distribution of the
%   selected minis-like events (all of them, unless TOP10PERCENTONLY is
%   set). MINIS is the summary output array produced by the DETECTMINIS
%   function and is described there. V is the electrophysiological
%   recording data vector in millivolts or nanoamperes. DT is the data
%   sampling interval in milliseconds (ms; scalar). TAU_M is the dendritic
%   passive membrane time constant (ms) estimated using impulses (scalar).
%   EXCLUDEDT is a vector containing excluded times (ms). TBLDURATION is a
%   scalar variable corresponding to the baseline duration in ms used for
%   estimating the amplitude of the average waveform. RTTYPE is a string
%   variable representing the rise time interval of choice: '10-90%' for a
%   10-90% rise time and '20-80%' for a 20-80% rise time. RISETIMEARRAYEXT
%   is a vector with the rise time range for classification extended by a
%   single element. DRAW is a logical variable with true corresponding to
%   presenting the graphical output of a function that includes the graph
%   of the average waveform. False corresponds to presenting no visual
%   output. WAVEFORM is a vector containing the average waveform. AVERAGEAMP
%   is a scalar variable corresponding to the mean amplitude of the
%   selected minis-like events. MEDIANAMP is a scalar variable
%   corresponding to the median amplitude of the selected minis-like
%   events. WAVEFORMRT is a rise time histogram vector of all detected
%   minis-like events (independent of TOP10PERCENTONLY, which only affects
%   which events are averaged into WAVEFORM/AVERAGEAMP/MEDIANAMP/
%   PARAMETERS.tau_m). PARAMETERS is a structure variable containing the
%   following fields:
%       'peak' - is a scalar corresponding to the peak value of the average
%           waveform in millivolts (mV) or nanoamperes (nA).
%       'BL' - is a scalar corresponding to the baseline value of the
%           average waveform, mV or nA.
%       'amplitude' - is a scalar corresponding to the amplitude value of
%           the average waveform, mV or nA.
%       'risetime' - is a scalar corresponding to the rise time value of
%           the average waveform, mV or nA.
%       'tau_m' - is a scalar corresponding to the passive membrane time
%           constant value of the average waveform in milliseconds.
%
%   [WAVEFORM WAVEFORMAMP WAVEFORMRT PARAMETERS F] = AVERAGEWAVEEFFECT(...)
%   In addition outputs F, the handle(s) of the average waveform graph and
%   (unless suppressed via DRAWRISETIMEFIGURE, see below) the rise-time
%   distribution graph.
%
%   AVERAGEWAVEEFFECT(..., draw, downGoing, top10PercentOnly,...
%       drawRiseTimeFigure) additionally accepts:
%       DOWNGOING (logical, default false) - minis(:,4) holds negative
%           amplitudes in the original trace convention, so the
%           largest-magnitude events and the waveform trough, not peak,
%           must be found.
%       TOP10PERCENTONLY (logical, default false) - restrict the average
%           to the top 10% largest-amplitude events instead of using all
%           of them. Used by the preprocessing pipeline (dataMinis.m).
%       DRAWRISETIMEFIGURE (logical, default true) - also draw the
%           second (rise-time-distribution) figure when DRAW is true. The
%           Detect/Detect-and-compare tasks pass false here, since that
%           histogram duplicates one already shown in their results panel;
%           WAVEFORMRT itself is still computed and returned regardless.
%


% Optional downGoing / top10PercentOnly / drawRiseTimeFigure flags:
if numel(varargin) >= 1 && ~isempty(varargin{1})
    downGoing = varargin{1};
else
    downGoing = false;
end
if numel(varargin) >= 2 && ~isempty(varargin{2})
    top10PercentOnly = varargin{2};
else
    top10PercentOnly = false;
end
if numel(varargin) >= 3 && ~isempty(varargin{3})
    drawRiseTimeFigure = varargin{3};
else
    drawRiseTimeFigure = true;
end

% Select which events go into the average waveform:
if top10PercentOnly
    % Restrict to the top 10% largest-amplitude events (used by the
    % preprocessing pipeline, dataMinis.m, for a cleaner, less
    % noise-contaminated waveform/tau_m estimate):
    if downGoing
        % Amplitudes are negative; pick the 10% with the most negative
        % values (i.e. the 10th percentile of amplitude = largest
        % magnitude):
        prctCutoff = prctile(minis(:,4), 10);
        iSelected = minis(:,4) <= prctCutoff;
    else
        prctCutoff = prctile(minis(:,4), 90);
        iSelected = minis(:,4) >= prctCutoff;
    end
else
    % Default: average all detected events:
    iSelected = true(size(minis,1), 1);
end
selectedPeakIndices = minis(iSelected,3);
averageAmp = mean(minis(iSelected,4));
medianAmp = median(minis(iSelected,4));

% Remove minis that overlap with excluded times:
iTau_m = round(tau_m/dt);
midWave = 8*iTau_m;
BLtime = 7*iTau_m;
minisToExclude = zeros(length(selectedPeakIndices), 1);
excludedInd = [1 round(excludedT/dt)];
for iMinis = 1:length(minisToExclude)
    distances = abs(excludedInd - selectedPeakIndices(iMinis));
    if min(distances) <= midWave
        minisToExclude(iMinis) = 1;
    end
end
selectedPeakIndices(logical(minisToExclude)) = [];
if isempty(selectedPeakIndices)
    waveform = [];
    waveformRT = [];
    parameters.tau_m = [];
    varargout = {[]};
    return
end

% Average selected minis:
iStart = selectedPeakIndices - midWave;
% Exclude any event whose extraction window would fall outside [1, length(V)].
% This can happen when large-amplitude events are very close to the start or
% end of the recording (a pre-existing edge case that is more easily reached
% after manual editing with the checkpoint resume feature).
iEndWave = selectedPeakIndices + midWave + 1;
validWindow = (iStart >= 1) & (iEndWave <= length(V));
if ~all(validWindow)
    nExcl = sum(~validWindow);
    disp(['averageWaveEffect: excluding ' num2str(nExcl) ' event(s) whose waveform ' ...
          'window exceeds the recording boundary.']);
    selectedPeakIndices = selectedPeakIndices(validWindow);
    iStart     = iStart(validWindow);
end
if isempty(selectedPeakIndices)
    waveform   = [];
    waveformRT = [];
    parameters.tau_m = [];
    varargout = {[]};
    return
end
waveformArray = repmat(iStart, 1, 2*midWave + 1);
waveformArray = waveformArray + repmat(1:2*midWave + 1, length(selectedPeakIndices), 1);
waveformArray = V(waveformArray);
waveform = mean(waveformArray,1);

% Estimate the peak, the baseline, the amplitude, and the 10-90% or 20-80% rise time:
if downGoing
    % Waveform is negative-going; find the minimum (trough) as the peak:
    [parameters.peak, iPeak] = min(waveform(midWave-100:midWave+1));
else
    [parameters.peak, iPeak] = max(waveform(midWave-100:midWave+1));
end
iPeak = midWave - 101 + iPeak;
halfBL = round(tBLduration/2/dt);
BL1 = waveform(iPeak-BLtime-halfBL+1:iPeak-BLtime+halfBL);
BL2 = waveform(iPeak+BLtime-halfBL:iPeak+BLtime+halfBL-1);
parameters.BL = (mean(BL1) + mean(BL2))/2;
parameters.amplitude = parameters.peak - parameters.BL;
waveTimeline = (1:length(waveform))*dt - dt;
if downGoing
    % calcRT always expects a positive-going waveform (peak > BL); pass
    % negated segment and negated baseline:
    [parameters.risetime, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, waveTimeline(1:iPeak), -waveform(1:iPeak), -parameters.BL);
else
    [parameters.risetime, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, waveTimeline(1:iPeak), waveform(1:iPeak), parameters.BL);
end

% Find the effective time constant:
waveformNorm = waveform - parameters.BL;
waveformNorm = gsmooth(waveformNorm, round(5/dt), 3);
if downGoing
    % For down-going events the normalised amplitude is negative; find the
    % point that exceeds 2/3 of the (negative) amplitude (i.e. the decay
    % has returned 1/3 of the way back toward BL):
    iStart = find(waveformNorm(iPeak:end) >= 2/3*parameters.amplitude, 1);
else
    iStart = find(waveformNorm(iPeak:end) <= 2/3*parameters.amplitude, 1);
end
if isempty(iStart)
    iStart = 1;
end
Vstart = waveformNorm(iPeak+iStart-1);
Vend = Vstart*1/exp(1);
if downGoing
    % Decay goes back toward 0 (less negative); find where it rises above Vend:
    iEnd = find(waveformNorm(iPeak:end) >= Vend, 1);
else
    iEnd = find(waveformNorm(iPeak:end) <= Vend, 1);
end
if isempty(iEnd)
    iEnd = length(waveformNorm)-iPeak+1;
end
parameters.tau_m = (iEnd - iStart)*dt;

% Smooth the data of interest:
waveSmooth = gsmooth(waveform, round(5/dt), 3);
Vstart = waveSmooth(iPeak+iStart-1);
Vend = waveSmooth(iPeak+iEnd-1);

% Draw a waveform graph:
BL = parameters.BL*ones(1,length(waveTimeline));
if draw
    f1 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83], 'Name', 'Averaged Minis-like Events', 'NumberTitle', 'off');
    plot(waveTimeline, waveform, 'k', 'markersize', 4);
    hold on
    xlim([waveTimeline(1) waveTimeline(end)]);
    if downGoing
        peakMarkWave = 'v';   % downward triangle for down-going events
    else
        peakMarkWave = '^';   % upward triangle for up-going events
    end
    plot(waveTimeline(iPeak), waveform(iPeak), peakMarkWave, 'markersize', 10, 'markeredgecolor', 'r');
    plot([t10 t50 t90], [waveform(i10) waveform(i50) waveform(i90)], 'r.', 'markersize', 10);
    plot(waveTimeline, BL, 'r');
    plot([waveTimeline(iPeak+iStart-1) waveTimeline(iPeak+iEnd-1)], [Vstart Vend], 'mx', 'markersize', 10);
    plot(waveTimeline(iPeak+iStart-1:iPeak+iEnd-1), waveSmooth(iPeak+iStart-1:iPeak+iEnd-1), 'm');
    title('Average waveform');
    xlabel('Time (ms)');
    ylabel('Voltage (mV)');
    if strcmpi(RTtype, '10-90%')
        RTlegend = '10,50,90% rise time';
    elseif strcmpi(RTtype, '20-80%')
        RTlegend = '20,50,80% rise time';
    end
    legend('Average waveform','Peak',RTlegend,'Baseline','Effective tau_m','Smoothed waveform','Location','NorthEast');
    hold off
    varargout = {f1};
end

% Rise time distribution of all detected events (independent of
% TOP10PERCENTONLY -- see the function header comment):
waveformRT = hist(minis(:,12), riseTimeArrayExt(2:end)); %#ok<*HIST>
waveformRT = [0 waveformRT(1:end-1)];
if draw && drawRiseTimeFigure
    if strcmpi(RTtype, '10-90%')
        figName = '10-90% Rise Time Distribution of Averaged Mini-like Events';
    elseif strcmpi(RTtype, '20-80%')
        figName = '20-80% Rise Time Distribution of Averaged Mini-like Events';
    end
    f2 = figure('Units', 'Normalized', 'Position', [0.005 0.06 0.99 0.83], 'Name', figName, 'NumberTitle', 'off');
    plot(riseTimeArrayExt(1:end-1), waveformRT,'b.-');
    hold on
    xlim([riseTimeArrayExt(1) riseTimeArrayExt(end-1)]);
    if strcmpi(RTtype, '10-90%')
        title('10-90% rise time distributions');
        xlabel('10-90% rise times (ms)');
    elseif strcmpi(RTtype, '20-80%')
        title('20-80% rise time distributions');
        xlabel('20-80% rise times (ms)');
    end
    ylabel('Number of Events');
    hold off
    varargout = {[f1 f2]};
end