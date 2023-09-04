function [waveform, parameters, varargout] = averageWaveform(minis, V, dt, SWduration, excludedT, tBLduration, peelDuration, cutoff, draw)
% AVERAGEWAVEFORM is a helper subfunction of detect minis. It estimates the
% average waveform of the minis detected using the detectMinis program.
%
%   [WAVEFORM PARAMETERS] = AVERAGEWAVEFORM(minis, V, dt, SWduration,...
%       excludedT, tBLduration, peelDuration, cutoff, draw)
%   Estimates the average waveform. MINIS is the matrix produced by the
%   minis detection algorithm of the detectMinis program. V is the
%   electrophysiological recording trace in milivolts (mV) or nanoamperes
%   (nA). DT is the sampling interval in miliseconds (ms). SWDURATION is
%   the beginning of the search window that corresponds to
%   searchParameters.SWstart as used by the minis detection algorithm, ms.
%   EXCLUDEDT is the vector of excluded times. TBLDURATION is the duration
%   of the baseline period, ms. PEELDURATION is a scalar variable
%   representing the duration of the period used for estimating the
%   membrane time constant, tau, in ms (or, in other words, for fitting the
%   regression line to the logarithm of the minis decay phase). CUTOFF is a
%   scalar variable representing the lower amplitude bound of minis used
%   for estimating the average waveform. Minis with smaller amplitudes than
%   whis value get rejected. DRAW is a string variable controling the
%   dislpay of the average waveform and regression line fitting fingures.
%   Setting 'on' turns on the figure displays, whereas setting to 'off'
%   turns them off.
%   WAVEFORM is a vector containg the average waveform, mV or nA.
%   PARAMETERS is a structure variable containing fields with average
%   waveform parameters 'peak' (mV or nA), 'BL' (baseline, mV or nA), 'Amp'
%   (amplitude, mV or nA), 'risetime' (ms), 'tau' (dendritic membrane time
%   constant, ms).
%
%   [WAVEFORM PARAMETERS F] = AVERAGEWAVEFORM(...)
%   In addition outputs the handles of the average waveform and regression
%   line fitting fingures. For this option to work, the user has to set the
%   input DRAW variable to 'on'.
%


minis(minis(:,4)<=cutoff,:) = [];
standardWaveformIndices = 0:12*round(SWduration/dt)-1;
standardWaveformIndices = repmat(standardWaveformIndices, length(minis(:,17)), 1);
waveformStartIndices = minis(:,17) - 3*SWduration/dt;

% Eliminate minis that contain excluded times:
midWave = ceil(size(standardWaveformIndices,2)/2);
waveformMidIndices = waveformStartIndices + midWave;
rowsToExclude = zeros(length(waveformMidIndices), 1);
excludedInd = [1 round(excludedT/dt)];
for iRow = 1:length(rowsToExclude)
    distances = abs(excludedInd - waveformMidIndices(iRow));
    if min(distances) <= midWave
        rowsToExclude(iRow) = 1;
    end
end
waveformStartIndices(logical(rowsToExclude)) = [];
clear waveformMidIndices rowsToExclude midWave excludedInd iRow distances minis SWduration excludedT

% Estimate the average waveform:
waveformArray = repmat(waveformStartIndices, 1, size(standardWaveformIndices, 2));
waveformArray = waveformArray + standardWaveformIndices(1:length(waveformStartIndices),:);
% waveformArray = V(waveformArray) - repmat(minis(logical(~rowsToExclude),20),1,size(waveformArray, 2)); % subtract the 'baseline'
waveformArray = V(waveformArray);
waveform = mean(waveformArray);
waveTimeline = 0:dt:length(waveform)*dt-dt;
clear waveformArray

% Estimate the peak, the baseline, and the amplitude:
[parameters.peak, iPeak] = max(waveform);
[~, iTrough] = min(waveform(1:iPeak));
iBLend = iTrough + round(.5*tBLduration/dt);
iBLstart = max([iTrough - round(.5*tBLduration/dt), 1]);
parameters.BL = mean(waveform(iBLstart:iBLend));
parameters.Amp = parameters.peak - parameters.BL;

% Calculate the rise-time:
[parameters.risetime, t10, i10, t50, i50, t90, i90] = calcRT(waveTimeline(1:iPeak), waveform(1:iPeak), parameters.BL);

% Visualise the waveform:
if strcmpi(draw,'on')
    f1 = figure('position', [50 50 1200 600], 'Name', 'Averaged Mini-like Events', 'NumberTitle', 'off');
    plot(waveTimeline, waveform, 'k', 'markersize', 4);
    hold on
    plot(waveTimeline(iPeak), waveform(iPeak), '^', 'markersize', 10, 'markeredgecolor', 'r');
    plot(t10, waveform(i10), 'r.', 'markersize', 10);
    plot(t50, waveform(i50), 'r.', 'markersize', 10);
    plot(t90, waveform(i90), 'r.', 'markersize', 10);
    plot(waveTimeline(iBLstart), parameters.BL, 'c.', 'markersize', 10);
    plot(waveTimeline(iBLend), parameters.BL, 'b.', 'markersize', 10);
    title('Smoothed Data');
    xlabel('Time (ms)');
    ylabel('Voltage (mV)');
    hold off
end


% Estimate the membrane time constant (tau):
subtrWaveform = waveform(iPeak:end) - parameters.BL;
iStartPeel = round(3/dt) + 1;
iEndPeel = iStartPeel-1 + round(peelDuration/dt);
lnPeel = real(log(subtrWaveform));
lnTimeline = (0:1:length(lnPeel)-1)*dt;
if exist('draw','var') && strcmpi(draw,'on')
    f2 = figure('position', [50 50 1200 600], 'Name', 'Natural Logarithm of the Decay Period', 'NumberTitle', 'off');
    plot(lnTimeline, lnPeel, 'k.', 'markersize', 6);
    hold on
    title('Smoothed Data');
    xlabel('Time (ms)');
    ylabel('Voltage (mV)');
    varargout = {[f1 f2]};
end

independent = [lnTimeline(iStartPeel:iEndPeel)' ones(iEndPeel-iStartPeel+1, 1)];
coefficients = regress(lnPeel(iStartPeel:iEndPeel)', independent);
regressionLine = coefficients(1)*lnTimeline(iStartPeel:iEndPeel) + coefficients(2);
if exist('draw','var') && strcmpi(draw,'on')
    plot(lnTimeline(iStartPeel:iEndPeel), regressionLine, 'r');
    hold off
end
parameters.tau = -1/coefficients(1);