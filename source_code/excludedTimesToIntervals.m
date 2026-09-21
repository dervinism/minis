function intervalsMs = excludedTimesToIntervals(excludedTimes)
% EXCLUDEDTIMESTOINTERVALS converts a targetExcludedTimes/noiseExcludedTimes
% -style struct (fields startPulse, endPulse, startGlitch, endGlitch, all
% in seconds -- the raw, GUI-level pulse/glitch definitions) into an Nx2
% matrix of [start end] excluded intervals in milliseconds, for use with
% INCIDENCERATE.
%
% Note: this uses the raw pulse/glitch boundaries, not the wider margin
% that detectMinis internally expands them by (to keep detection clear of
% search-window/baseline artefacts near each pulse). The difference is
% typically tens of milliseconds against a 1-second bin, so it is ignored
% here for simplicity.

starts = [];
ends = [];
if isfield(excludedTimes, 'startPulse') && isfield(excludedTimes, 'endPulse')
    starts = [starts; excludedTimes.startPulse(:)];
    ends = [ends; excludedTimes.endPulse(:)];
end
if isfield(excludedTimes, 'startGlitch') && isfield(excludedTimes, 'endGlitch')
    starts = [starts; excludedTimes.startGlitch(:)];
    ends = [ends; excludedTimes.endGlitch(:)];
end

valid = ~isnan(starts) & ~isnan(ends);
intervalsMs = 1000 * [starts(valid) ends(valid)];
end
