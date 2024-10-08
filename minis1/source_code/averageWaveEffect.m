function [waveform, averageAmp, medianAmp, waveformRT, parameters, varargout] = averageWaveEffect(minis, V, dt, tau_m, excludedT,...
  tBLduration, RTtype, riseTimeArrayExt, draw)
% AVERAGEWAVEEFFECT is a helper function of DETECTMINIS. It estimates the
% average waveform of the 10% largest minis-like events as detected using
% the MINIS program.
%
%   [WAVEFORM WAVEFORMAMP WAVEFORMRT PARAMETERS] = AVERAGEWAVEEFFECT(...
%       minis, V, dt, tau_m, excludedT, tBLduration, RTtype,...
%       riseTimeArray, draw)
%   Estimates the average waveform and the rise time distribution of the
%   10% largest minis-like events. MINIS is the summary output array
%   produced by the DETECTMINIS function and is described there. V is the
%   electrophysiological recording data vector in millivolts or nanoamperes.
%   DT is the data sampling interval in milliseconds (ms; scalar). TAU_M is
%   the dendritic passive membrane time constant (ms) estimated using
%   impulses (scalar). EXCLUDEDT is a vector containing excluded times (ms).
%   TBLDURATION is a scalar variable corresponding to the baseline duration
%   in ms used for estimating the amplitude of the average waveform. RTTYPE
%   is a string variable representing the rise time interval of choice:
%   '10-90%' for a 10-90% rise time and '20-80%' for a 20-80% rise time.
%   RISETIMEARRAYEXT is a vector with the rise time range for
%   classification extended by a single element. DRAW is a logical variable
%   with true corresponding to presenting the graphical output of a
%   function that includes the graphs of the average waveform and the rise
%   time distribution. False corresponds to presenting no visual output.
%   WAVEFORM is a vector containing the average waveform. AVERAGEAMP is a
%   scalar variable corresponding to the mean amplitude of the top 10% of
%   the detected largest minis-like events. MEDIANAMP is a scalar variable
%   corresponding to the median amplitude of the top 10% of the detected
%   largest minis-like events. WAVEFORMRT is a rise time histogram vector
%   of the top 10% of the detected largest minis-like events. PARAMETERS is
%   a structure variable containing the following fields:
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
%   In addition outputs handles corresponding to the average waveform graph
%   and the rise time histogram of the top 10% largest minis-like events.
%


% Find minis larger than the 90th percentile of minis amplitudes:
prct90 = prctile(minis(:,4),90);
iPrct90 = minis(:,4) >= prct90;
minis90prc = minis(iPrct90,3);
averageAmp = mean(minis(iPrct90,4));
medianAmp = median(minis(iPrct90,4));

% Remove minis that overlap with excluded times:
iTau_m = round(tau_m/dt);
midWave = 8*iTau_m;
BLtime = 7*iTau_m;
minisToExclude = zeros(length(minis90prc), 1);
excludedInd = [1 round(excludedT/dt)];
for iMinis = 1:length(minisToExclude)
    distances = abs(excludedInd - minis90prc(iMinis));
    if min(distances) <= midWave
        minisToExclude(iMinis) = 1;
    end
end
minis90prc(logical(minisToExclude)) = [];
if isempty(minis90prc)
    waveform = [];
    waveformRT = [];
    parameters.tau_m = [];
    varargout = {[]};
    return
end

% Average selected minis:
iStart = minis90prc - midWave;
waveformArray = repmat(iStart, 1, 2*midWave + 1);
waveformArray = waveformArray + repmat(1:2*midWave + 1, length(minis90prc), 1);
waveformArray = V(waveformArray);
waveform = mean(waveformArray,1);

% Estimate the peak, the baseline, the amplitude, and the 10-90% or 20-80% rise time:
[parameters.peak, iPeak] = max(waveform(midWave-100:midWave+1));
iPeak = midWave - 101 + iPeak;
halfBL = round(tBLduration/2/dt);
BL1 = waveform(iPeak-BLtime-halfBL+1:iPeak-BLtime+halfBL);
BL2 = waveform(iPeak+BLtime-halfBL:iPeak+BLtime+halfBL-1);
parameters.BL = (mean(BL1) + mean(BL2))/2;
parameters.amplitude = parameters.peak - parameters.BL;
waveTimeline = (1:length(waveform))*dt - dt;
[parameters.risetime, t10, i10, t50, i50, t90, i90] = calcRT(RTtype, waveTimeline(1:iPeak), waveform(1:iPeak), parameters.BL);

% Find the effective time constant:
waveformNorm = waveform - parameters.BL;
waveformNorm = gsmooth(waveformNorm, round(5/dt), 3);
iStart = find(waveformNorm(iPeak:end) <= 2/3*parameters.amplitude, 1);
if isempty(iStart)
    iStart = 1;
end
Vstart = waveformNorm(iPeak+iStart-1);
Vend = Vstart*1/exp(1);
iEnd = find(waveformNorm(iPeak:end) <= Vend, 1);
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
    plot(waveTimeline(iPeak), waveform(iPeak), '^', 'markersize', 10, 'markeredgecolor', 'r');
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
end

% Calculate the 10-90% or 20-80% rise time distribution:
waveformRT = hist(minis(:,12), riseTimeArrayExt(2:end)); %#ok<*HIST>
waveformRT = [0 waveformRT(1:end-1)];
if draw
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
        title('10-90% Rise Time Distribution of Averaged Mini-like Events');
        title('10-90% rise time distributions');
        xlabel('10-90% rise times (ms)');
    elseif strcmpi(RTtype, '20-80%')
        title('20-80% Rise Time Distribution of Averaged Mini-like Events');
        title('20-80% rise time distributions');
        xlabel('20-80% rise times (ms)');
    end
    ylabel('Number of Events');
    hold off
    varargout = {[f1 f2]};
end