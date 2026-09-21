function rates = incidenceRate(peakTimesMs, recordingDurationMs, excludedIntervalsMs)
% INCIDENCERATE computes, for each detected event, the incidence rate
% (events per second) of the 1-second window containing that event's
% peak, assigned retrospectively to every event within that window.
%
%   RATES = INCIDENCERATE(peakTimesMs, recordingDurationMs, excludedIntervalsMs)
%   bins events into consecutive 1-second windows starting at t=0 (i.e.
%   [0,1), [1,2), ...). Each window's rate is:
%       (number of events whose peak falls in the window)
%       / (effective duration of the window, in seconds)
%   The effective duration accounts for:
%     - the recording's actual length: the final window is normally
%       shorter than 1 second and is divided by its true duration rather
%       than assumed to be a full second;
%     - excluded time regions (pulses/glitches) overlapping the window,
%       which are subtracted from its nominal duration.
%   A window whose effective duration is zero (entirely excluded, or
%   beyond the recording) yields NaN for any event assigned to it -- this
%   should not normally occur since events are not detected in excluded
%   regions, but is handled defensively.
%
%   PEAKTIMESMS is a vector of detected-event peak times (ms).
%   RECORDINGDURATIONMS is the total recording length (ms, scalar).
%   EXCLUDEDINTERVALSMS is an Nx2 matrix of [start end] excluded time
%   intervals (ms); pass an empty array ([] or zeros(0,2)) if none.
%
%   Each element of RATES corresponds to the same-indexed element of
%   PEAKTIMESMS.

binWidthMs = 1000;
peakTimesMs = peakTimesMs(:);

if isempty(peakTimesMs)
    rates = zeros(0,1);
    return
end

if nargin < 3 || isempty(excludedIntervalsMs)
    excludedIntervalsMs = zeros(0,2);
end

nBins = max(1, ceil(recordingDurationMs / binWidthMs));

% Nominal duration of each bin (s); the last bin is clipped to the
% recording's actual end and so may be shorter than 1 s:
binStartsMs = (0:nBins-1)' * binWidthMs;
binEndsMs = min(binStartsMs + binWidthMs, recordingDurationMs);
nominalDurationS = (binEndsMs - binStartsMs) / 1000;

% Excluded duration overlapping each bin (s), via interval overlap:
excludedDurationS = zeros(nBins, 1);
for iInt = 1:size(excludedIntervalsMs, 1)
    overlapStart = max(binStartsMs, excludedIntervalsMs(iInt,1));
    overlapEnd = min(binEndsMs, excludedIntervalsMs(iInt,2));
    overlapMs = max(overlapEnd - overlapStart, 0);
    excludedDurationS = excludedDurationS + overlapMs / 1000;
end

effectiveDurationS = max(nominalDurationS - excludedDurationS, 0);

% Event count per bin:
eventBinIdx = floor(peakTimesMs / binWidthMs) + 1;
eventBinIdx = min(max(eventBinIdx, 1), nBins);   % clamp rounding edge cases
eventCounts = accumarray(eventBinIdx, 1, [nBins, 1]);

binRates = nan(nBins, 1);
hasData = effectiveDurationS > 0;
binRates(hasData) = eventCounts(hasData) ./ effectiveDurationS(hasData);

rates = binRates(eventBinIdx);
end
